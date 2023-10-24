## Bart
#' @useDynLib cppspBART1
#' @importFrom Rcpp sourceCpp
#'

# Getting the BART wrapped function
#' @export
spbart <- function(x_train,
                  y,
                  x_test,
                  n_tree = 2,
                  node_min_size = 5,
                  n_mcmc = 2000,
                  n_burn = 500,
                  alpha = 0.95,
                  beta = 2,
                  df = 3,
                  sigquant = 0.9,
                  kappa = 2,
                  tau = 100,
                  scale_bool = TRUE,
                  stump = FALSE,
                  no_rotation_bool = FALSE,
                  numcut = 100L
                  ) {


        # Interaction term verification
        if(interaction_term & is.null(interaction_list)){
                stop("Define the interaction list.")
        }

        # Verifying if x_train and x_test are matrices
        if(!is.data.frame(x_train) || !is.data.frame(x_test)){
                stop("Insert valid data.frame for both data and xnew.")
        }


        # Getting the valid
        dummy_x <- base_dummyVars(x_train)

        # Create a list
        if(length(dummy_x$facVars)!=0){
                for(i in 1:length(dummy_x$facVars)){
                        # See if the levels of the test and train matches
                        if(!all(levels(x_train[[dummy_x$facVars[i]]])==levels(x_test[[dummy_x$facVars[i]]]))){
                                levels(x_test[[dummy_x$facVars[[i]]]]) <- levels(x_train[[dummy_x$facVars[[i]]]])
                        }
                        df_aux <- data.frame( x = x_train[,dummy_x$facVars[i]],y)
                        formula_aux <- stats::aggregate(y~x,df_aux,mean)
                        formula_aux$y <- rank(formula_aux$y)
                        x_train[[dummy_x$facVars[i]]] <- as.numeric(factor(x_train[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

                        # Doing the same for the test set
                        x_test[[dummy_x$facVars[i]]] <- as.numeric(factor(x_test[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

                }
        }


        # Getting the train and test set
        x_train_scale <- as.matrix(x_train)
        x_test_scale <- as.matrix(x_test)

        # Scaling x
        x_min <- apply(as.matrix(x_train_scale),2,min)
        x_max <- apply(as.matrix(x_train_scale),2,max)

        # Storing the original
        x_train_original <- x_train
        x_test_original <- x_test


        # Normalising all the columns
        for(i in 1:ncol(x_train)){
                x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
                x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
        }



        # Creating the numcuts matrix of splitting rules
        xcut_m <- matrix(NA,nrow = numcut,ncol = ncol(x_train_scale))
        for(i in 1:ncol(x_train_scale)){

                if(nrow(x_train_scale)<numcut){
                        xcut_m[,i] <- sort(x_train_scale[,i])
                } else {
                        xcut_m[,i] <- seq(min(x_train_scale[,i]),
                                          max(x_train_scale[,i]),
                                          length.out = numcut+2)[-c(1,numcut+2)]
                }
        }




        # =========================================================================================================
        # Getting the Splines Basis functions
        # =========================================================================================================
        # Creating a list of basis functions
        B_train_obj <- B_test_obj <- vector("list",length = length(dummy_x$continuousVars))

        # Setting new parameters for the spline
        ndx <- nIknots+1
        ord_ <- 4
        degree_ <- 3
        x_min_sp <- apply(x_train_scale,2,min)
        x_max_sp <- apply(x_train_scale,2,max)
        dx <- (x_max_sp-x_min_sp)/ndx

        # New_knots
        new_knots <- matrix()
        new_knots <- matrix(mapply(x_min_sp,x_max_sp,dx, FUN = function(MIN,MAX,DX){seq(from = MIN-(ord_-1)*DX, to = MAX+(ord_-1)*DX, by = DX)}), ncol = length(dummy_x$continuousVars)) # MIN and MAX are 0 and 1 respectively, because of the scale
        colnames(new_knots) <- dummy_x$continuousVars

        # Selecting which one gonna be used bs or splines.des()
        if(interaction_term){

                D_train <- matrix(NA,
                                  nrow = nrow(x_train_scale),
                                  ncol = (nrow(new_knots)-ord_)*(length(dummy_x$continuousVars))+(nrow(new_knots)-ord_)^2*length(interaction_list))

                D_test <- matrix(NA,
                                 nrow = nrow(x_test_scale),
                                 ncol = (nrow(new_knots)-ord_)*(length(dummy_x$continuousVars))+(nrow(new_knots)-ord_)^2*length(interaction_list))
        } else {
                D_train <- matrix(NA,
                                  nrow = nrow(x_train_scale),
                                  ncol = (nrow(new_knots)-ord_)*length(dummy_x$continuousVars))

                D_test <- matrix(NA,
                                 nrow = nrow(x_test_scale),
                                 ncol = (nrow(new_knots)-ord_)*length(dummy_x$continuousVars))
        }


        # Selecting the basis size.
        if(interaction_term){ # Checking the interaction
                basis_size <- (nrow(new_knots)-ord_)     # Change this value to the desired size of each sublist
        } else {
                basis_size <- (nrow(new_knots)-ord_)     # Change this value to the desired size of each sublist
        }


        D_seq <- 1:ncol(D_train)  # Replace this with the columns of D

        # Creating a vector
        basis_subindex <- split(D_seq, rep(1:((length(D_seq) %/% basis_size)), each = basis_size, length.out = length(D_seq)))

        start_ <- NCOL(x_train_scale)+1

        # Iterate over pairs from the start_index
        for (i in 1:(length(interaction_list))) {
                basis_subindex[[start_]] <- unlist(basis_subindex[start_:length(basis_subindex)])
                basis_subindex[(start_+1):length(basis_subindex)] <- NULL
        }

        # Creating the natural B-spline for each predictor
        for(i in 1:NCOL(x_train_scale)){


                # Modify the basis only with respect to the main effects at the moment
                B_train_obj[[i]] <- splines::spline.des(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                                        knots = new_knots[,dummy_x$continuousVars[i]],
                                                        ord = ord_,
                                                        derivs = 0*x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design


                # Adding interaction

                # Returning to MOTR-BART
                if(length(basis_subindex[[i]])!= ncol(B_train_obj[[i]])){
                        stop("Error on the basis generation")
                }

                D_train[,basis_subindex[[i]]] <- as.matrix(B_train_obj[[i]])


                # For the test setting

                B_test_obj[[i]] <- splines::spline.des(x = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                                       knots = new_knots[,dummy_x$continuousVars[i]],
                                                       ord = ord_,
                                                       derivs = 0*x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design


                # Returning to MOTR-BART
                if(length(basis_subindex[[i]])!= ncol(B_test_obj[[i]])){
                        stop("Error on the basis generation")
                }

                D_test[,basis_subindex[[i]]] <- as.matrix(B_test_obj[[i]])

        }

        # Interaction matrix list to be used in the penalised
        interaction_matrix_list <- vector("list",length = length(interaction_list))

        # Adding the interaction basis
        if(interaction_term){
                for (jj in 1:length(interaction_list)) {
                        interaction_matrix_list[[jj]] <- multiply_matrices_general(A = B_train_obj[[interaction_list[[jj]][1]]],B = B_train_obj[[interaction_list[[jj]][2]]])
                        D_train[,basis_subindex[[NCOL(x_train_scale)+jj]]] <- interaction_matrix_list[[jj]]
                        D_test[,basis_subindex[[NCOL(x_test_scale)+jj]]] <- multiply_matrices_general(A = B_test_obj[[interaction_list[[jj]][1]]],B = B_test_obj[[interaction_list[[jj]][2]]])
                }
        }

        # ==== COMMENTTED FUNCTIONS BELOW NOT RUN WITH IF NOT INSIDE THE FUNCTION
        # selected_var <- 1
        # D_subset <- D_train[,basis_subindex[[selected_var]]]
        # plot(NULL,ylim = range(D_subset), xlim = range(x_train_scale[,selected_var]), main = "Use BS: FALSE")
        # for(i in 1:ncol(D_subset)){
        #   points(x_train_scale[,selected_var], D_subset[,i], pch = 20, col = ggplot2::alpha(i,0.5))
        # }

        # Adapting to go the MOTR-BART appraoch.
        if(motrbart_bool){
                D_train <- x_train_scale
                D_test <- x_test_scale

                basis_size <- 1 # Change this value to the desired size of each sublist
                D_seq <- 1:ncol(D_train)  # Replace this with the columns of D

                # Creating a vector
                basis_subindex <- split(D_seq, rep(1:(length(D_seq) %/% basis_size), each = basis_size, length.out = length(D_seq)))
        }

        # Scaling the y
        min_y <- min(y_train)
        max_y <- max(y_train)

        # Getting the min and max for each column
        min_x <- apply(x_train_scale,2,min)
        max_x <- apply(x_train_scale, 2, max)

        # Scaling "y"
        if(scale_bool){
                y_scale <- normalize_bart(y = y_train,a = min_y,b = max_y)

                # New update
                m_tilda <- mean(diag(tcrossprod(D_train)))
                # Maybe need to change that in the future
                tau_mu <- 4*n_tree*(kappa^2)*(m_tilda)

        } else {
                y_scale <- y_train
                # New parameter update
                m_tilda <- mean(diag(tcrossprod(D_train)))
                tau_mu <- (4*n_tree*(kappa^2)*(m_tilda))/((max_y-min_y)^2)
        }


        # Getting the naive sigma value
        nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

        # Calculating tau hyperparam
        a_tau <- df/2

        # Calculating lambda
        qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
        lambda <- (nsigma*nsigma*qchi)/df
        d_tau <- (lambda*df)/2


        # Getting hyperparameters for \tau_beta_j
        a_tau_beta_j <- df/2
        sigquant_beta <- 0.99
        nsigma_beta <- tau_mu^(-1/2)

        # Calculating lambda
        qchi_beta <- stats::qchisq(p = 1-sigquant_beta,df = df,lower.tail = 1,ncp = 0)
        lambda_beta <- (nsigma_beta*nsigma_beta*qchi_beta)/df
        d_tau_beta_j <- (lambda_beta*df)/2


        # Call the bart function
        tau_init <- nsigma^(-2)
        mu_init <- mean(y_scale)


        # Creating a list to store all the main effects for the sum of trees
        if(main_effects_pred){

                # Calculating the main effects
                if(interaction_term){

                        # Interacting
                        main_effects_train_list <- main_effects_test_list <- vector("list", length = length(dummy_x$continuousVars)+length(interaction_list))

                        for(list_size in 1:length(main_effects_train_list)){
                                main_effects_train_list[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_train))
                                main_effects_test_list[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_test))
                        }

                        # Renaming the main list
                        if(interaction_term){
                                internames_ <- lapply(interaction_list,function(vars_){paste0("x",paste0(vars_,collapse = ":"))})
                                names(main_effects_train_list) <- names(main_effects_test_list) <- c(dummy_x$continuousVars, internames_)
                        } else {
                                names(main_effects_train_list) <- names(main_effects_test_list) <-dummy_x$continuousVars
                        }


                } else { # In case there are no interactions
                        main_effects_train_list <- main_effects_test_list <- vector("list", length = length(dummy_x$continuousVars))

                        for(list_size in 1:length(dummy_x$continuousVars)){
                                main_effects_train_list[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_train))
                                main_effects_test_list[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_test))
                        }

                        names(main_effects_train_list) <- names(main_effects_test_list) <- dummy_x$continuousVars
                }
        } else {
                main_effects_train_list <- main_effects_test_list <- NULL
        }

        # Gonna create a list of lists to store all the indexes for all split rules and cutpoints
        all_var_splits <- vector("list",ncol(x_train_scale))
        names(all_var_splits) <- colnames(x_train_scale)



        # Iterating over all possible x.columns
        for(i in 1:length(all_var_splits)){

                # Creating the dummy for a list of index to store all numeric split values
                all_cut_points <- vector("list", nrow(xcut_m))


                for(j in 1:length(all_cut_points)){

                        # Getting the node indexes object
                        left_train_list <- vector("list",length = 1L)
                        names(left_train_list) <- "left_train"
                        right_train_list <- vector("list",length = 1L)
                        names(right_train_list) <- "right_train"
                        left_test_list <- vector("list",length = 1L)
                        names(left_test_list) <- "left_test"
                        right_test_list <- vector("list",length = 1L)
                        names(right_test_list) <- "right_test"

                        node_index <- append(left_train_list, right_train_list) |>
                                append(left_test_list) |> append(right_test_list)

                        all_cut_points[[j]]$left_train <-  which(x_train_scale[,i] < xcut_m[j,i])
                        all_cut_points[[j]]$right_train <-  which(x_train_scale[,i] >= xcut_m[j,i])
                        all_cut_points[[j]]$left_test <-  which(x_test_scale[,i] < xcut_m[j,i])
                        all_cut_points[[j]]$right_test <-  which(x_test_scale[,i] >= xcut_m[j,i])

                }

                all_var_splits[[i]] <- all_cut_points

        }

        # Creating the penalty matrix

        all_P <- replicate(NCOL(x_train_scale),
                           P_gen(D_train_ = B_train_obj[[1]],dif_order_ = dif_order,tau_mu_ = 1),simplify = FALSE)
        if(interaction_term){
                # Adding the penalty term for the interactions
                for( ii in 1:length(interaction_list)){
                        all_P_aux <- list(kronecker(all_P[[interaction_list[[ii]][1]]],all_P[[interaction_list[[ii]][2]]]))
                }

                all_P <- append(all_P,all_P_aux)
        }

        P_train <- as.matrix(Matrix::bdiag(all_P))

        # Generating the BART obj
        bart_obj <- cppbart(x_train_scale,
                                y_scale,
                                x_test_scale,
                                xcut_m,
                                n_tree,
                                node_min_size,
                                n_mcmc,
                                n_burn,
                                tau_init,
                                mu_init,
                                tau_mu,
                                alpha,
                                beta,
                                a_tau,d_tau,
                                stump)



        # ==== Scaling back the data

        if(scale_bool){
             # Tidying up the posterior elements
             y_train_post <- unnormalize_bart(z = bart_obj[[1]],a = min_y,b = max_y)
             y_test_post <- unnormalize_bart(z = bart_obj[[2]],a = min_y,b = max_y)
             for(i in 1:round(n_mcmc-n_burn)){
                     all_tree_post[[i]] <-  unnormalize_bart(z = bart_obj[[4]][,,i],a = min_y,b = max_y)
             }
             tau_post <- bart_obj[[3]]/((max_y-min_y)^2)
             all_tau_post <- bart_obj[[7]]/((max_y-min_y)^2)
        } else {
             y_train_post <- bart_obj[[1]]
             y_test_post <- bart_obj[[2]]
             tau_post <- bart_obj[[3]]
             for(i in 1:round(n_mcmc-n_burn)){
                     all_tree_post[[i]] <-  bart_obj[[4]][,,i]
             }
             all_tau_post <- bart_obj[[7]]


        }

     # Return the list with all objects and parameters
     return(list(y_hat = y_train_post,
                 y_hat_test = y_test_post,
                 tau_post = tau_post,
                 all_tau_post = all_tau_post,
                 all_tree_post = all_tree_post,
                 prior = list(n_tree = n_tree,
                              alpha = alpha,
                              beta = beta,
                              tau_mu = tau_mu,
                              a_tau = a_tau,
                              d_tau = d_tau),
                 mcmc = list(n_mcmc = n_mcmc,
                             n_burn = n_burn),
                 data = list(x_train = x_train,
                             y = y,
                             x_test = x_test,
                             move_proposal = bart_obj[[5]],
                             move_acceptance = bart_obj[[6]])))
}


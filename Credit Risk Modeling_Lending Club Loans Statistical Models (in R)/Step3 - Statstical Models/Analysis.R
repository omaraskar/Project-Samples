########################### Analysis file ###################
# Input for this code is the tidy data set "LC_data"

# Load the tidy data (should be saved as an R file "LC_data", adjust path as needed)
load("c:/users/david/desktop/lending club data/LC_data")

# Computer with 16 GB ram - may be needed as many steps are computationally intensive
# R version 3.3.2 (10-31-2016 release)
# Load required libraries, install packages as needed
library(dplyr, warn.conflicts = FALSE) # v 0.5.0
library(tidyr, warn.conflicts = FALSE) # v 0.6.0
library(corrplot, warn.conflicts = FALSE) # v 0.77, for visualizing correlation matrices
library(MASS, warn.conflicts = FALSE) # v 7.3-45, for discriminant analysis

# Remove location variables,employment titles, ids, and dates
# These variables will not be used in the predicitve modeling
# Also remove variables that were used to create
# adjusted income, as these are linear combinations
# of the original income variables and will cause issues in modeling (multicolinearity)
LC_data = LC_data[,-(which(names(LC_data) %in% c("emp_title","zip_code","addr_state","grade",
                 "state_zip_code","county_id","county_name",
                 "county_latitude","county_longitude","population_2000",
                 "id","member_id","issue_d","earliest_cr_line",
                 "loan_status","annual_inc","pct_loan_income","median_rent_2_bed",
                 "loan_amnt", "pub_rec")))]

# Visualize the correlation matrix to look for 
# possible multicolinearity
M = cor(LC_data[,-c(4,6)], use = "complete.obs")
corrplot(M, method = "ellipse", tl.col = "black", tl.cex = 0.8)

# Roughly 25,000 observations of the variables
# bc_util - percent_bc_gt_75 are missing 
# From the correlation plot, these variables were 
# we highly correlated with each other and also
# with revol_util, inq_last_6mths, open_acc
# These variables will be dropped for the logistic regression model building, as
# most of the information they carry is captured by revol_util,inq_last_6mths, and open_acc

# Remove these variables from the data
LC_complete = LC_data[,-(16:26)]

# Remove missing data
LC_complete = LC_complete[complete.cases(LC_complete),]

# Visualize correlations of numeric variables in complete data set
C = cor(LC_complete[,-c(4,6)])
corrplot(C, method = "ellipse", type = "lower",tl.col = "black", tl.cex = 0.8)

# Obtain summary statistics for the null model
Null_model_stats = LC_complete %>% group_by(loan_risk) %>% 
  summarise(Non_default = sum(I(loan_default == 0)),
            Default = sum(I(loan_default == 1))) %>% 
  ungroup() %>% 
  mutate(Total_loans = Non_default + Default,
         Default_rate = round(100*(Default/Total_loans),1))

################ Logistic Regression ####################################

# Perform variable selection using stepwise logistic regression
# Logistic Model with all variables
Full_model = glm(loan_default ~ ., data = LC_complete, family = "binomial")

#Logistic model with only intercept term
Null_model = glm(loan_default ~ 1, data = LC_complete, family = "binomial")

# Perform mixed stepwise regression to select best variables
# Set seed for reproducable results
set.seed(243434)
Selected_model = step(Null_model, scope = list(lower = formula(Null_model),
                      upper = formula(Full_model)),
                      direction = "both", trace = 0)

############## Cross Validation for Logistic Model ####################
# Create training and test data sets
# Test model on 1000 random test data sets

# Set seed for reproducable results
set.seed(1456)

# Set up matrices and vectors for data collection 
training_index = numeric()
test_index = numeric()
misclass_rates = matrix(NA, nrow = 1000, ncol = 4)
correct_default_rates = matrix(NA, nrow = 1000, ncol = 4)
expected_error_rates = matrix(NA, nrow = 1000, ncol = 4)

# Test the model on 1000 random test sets and store the error and true positive rates
for (i in 1:1000) {
training_index = sample(nrow(LC_complete),65000)
test_index = -training_index

LC_train = LC_complete[training_index,]
LC_test = LC_complete[test_index,]

# Obtain Error rates for best logistic model
LC_test$Logistic_probs = predict(Selected_model, LC_test[-1],type = "response")
LC_test$Logistic_predict = ifelse(LC_test$Logistic_probs >= 0.20, 1, 0)

# Confusion Matrices by Group
All = with(LC_test,table(loan_default,Logistic_predict))
Low = table(LC_test[LC_test$loan_risk=="Low","loan_default"],
            LC_test[LC_test$loan_risk=="Low","Logistic_predict"])
Moderate = table(LC_test[LC_test$loan_risk=="Moderate","loan_default"],
            LC_test[LC_test$loan_risk=="Moderate","Logistic_predict"])
High = table(LC_test[LC_test$loan_risk=="High","loan_default"],
            LC_test[LC_test$loan_risk=="High","Logistic_predict"])

# Misclassification rate by risk
misclass_rates[i,] = c(round(100-(100*(sum(diag(All))/sum(All))),1),
                       round(100-(100*(sum(diag(Low))/sum(Low))),1),
                       round(100-(100*(sum(diag(Moderate))/sum(Moderate))),1),
                       round(100-(100*(sum(diag(High))/sum(High))),1))

# True Defaults Correctly Identified
correct_default_rates[i,] = c(round((100*(All[2,2]/sum(All[2,]))),1),
                               round((100*(Low[2,2]/sum(Low[2,]))),1),
                               round((100*(Moderate[2,2]/sum(Moderate[2,]))),1),
                               round((100*(High[2,2]/sum(High[2,]))),1))

# Store expected error rates (# bad loans/ # predicited good loans)
expected_error_rates[i,] = c(round((100*(All[2,1]/sum(All[,1]))),1),
                              round((100*(Low[2,1]/sum(Low[,1]))),1),
                              round((100*(Moderate[2,1]/sum(Moderate[,1]))),1),
                              round((100*(High[2,1]/sum(High[,1]))),1))
}

# Save output from cross validation in data frames for graphics
misclass_rates_logistic = data.frame(misclass_rates)
names(misclass_rates_logistic) = c("All_loans","Low_risk","Moderate_risk","High_risk")

correct_default_rates_logistic = data.frame(correct_default_rates)
names(correct_default_rates_logistic) = c("All_loans","Low_risk","Moderate_risk","High_risk")

expected_error_rates_logistic = data.frame(expected_error_rates)
names(expected_error_rates_logistic) = c("All_loans","Low_risk","Moderate_risk","High_risk")

# Transpose results for tidy data
misclass_rates_logistic = gather(misclass_rates_logistic,
                                 key = "Loan_Risk",
                                 value = "Estimated_Misclassification_Rate",
                                 1:4)
misclass_rates_logistic$Loan_Risk = 
  factor(misclass_rates_logistic$Loan_Risk,
         levels = c("All_loans","Low_risk","Moderate_risk","High_risk"),
         labels = c("All_loans","Low_risk","Moderate_risk","High_risk"))

correct_default_rates_logistic = gather(correct_default_rates_logistic,
                                 key = "Loan_Risk",
                                 value = "Estimated_Identitfy_True_Default_Rate",
                                 1:4)
correct_default_rates_logistic$Loan_Risk = 
  factor(correct_default_rates_logistic$Loan_Risk,
         levels = c("All_loans","Low_risk","Moderate_risk","High_risk"),
         labels = c("All_loans","Low_risk","Moderate_risk","High_risk"))

expected_error_rates_logistic = gather(expected_error_rates_logistic,
                                 key = "Loan_Risk",
                                 value = "Estimated_Expected_Error_Rate",
                                 1:4)
expected_error_rates_logistic$Loan_Risk = 
  factor(expected_error_rates_logistic$Loan_Risk,
         levels = c("All_loans","Low_risk","Moderate_risk","High_risk"),
         labels = c("All_loans","Low_risk","Moderate_risk","High_risk"))


# Save the results for later use
save(misclass_rates_logistic, 
     file = "c:/users/dsvancer/desktop/r directory/lending club data/misclass_rates_logistic")

save(correct_default_rates_logistic, 
     file = "c:/users/dsvancer/desktop/r directory/lending club data/correct_default_rates_logistic")

save(expected_error_rates_logistic, 
     file = "c:/users/dsvancer/desktop/r directory/lending club data/expected_error_rates_logistic")


# Summarize results
Misclassification_results = misclass_rates_logistic %>% 
                    group_by(Loan_Risk) %>% 
                    summarise(Min = min(Estimated_Misclassification_Rate),
                    p25 = quantile(Estimated_Misclassification_Rate,probs = 0.25),
                    Median = median(Estimated_Misclassification_Rate),
                    p75 = quantile(Estimated_Misclassification_Rate,probs = 0.75),
                    Max = max(Estimated_Misclassification_Rate),
                    Mean = mean(Estimated_Misclassification_Rate),
                    SD = sd(Estimated_Misclassification_Rate))
                            
Correct_default_results = correct_default_rates_logistic %>% group_by(Loan_Risk) %>% 
  summarise(Min = min(Estimated_Identitfy_True_Default_Rate),
            p25 = quantile(Estimated_Identitfy_True_Default_Rate,
                           probs = 0.25),
            Median = median(Estimated_Identitfy_True_Default_Rate),
            p75 = quantile(Estimated_Identitfy_True_Default_Rate,
                           probs = 0.75),
            Max = max(Estimated_Identitfy_True_Default_Rate),
            Mean = mean(Estimated_Identitfy_True_Default_Rate),
            SD = sd(Estimated_Identitfy_True_Default_Rate))

Error_rate_results = expected_error_rates_logistic %>% group_by(Loan_Risk) %>% 
  summarise(Min = min(Estimated_Expected_Error_Rate),
            p25 = quantile(Estimated_Expected_Error_Rate,
                           probs = 0.25),
            Median = median(Estimated_Expected_Error_Rate),
            p75 = quantile(Estimated_Expected_Error_Rate,
                           probs = 0.75),
            Max = max(Estimated_Expected_Error_Rate),
            Mean = mean(Estimated_Expected_Error_Rate),
            SD = sd(Estimated_Expected_Error_Rate))


# Visualize CV results
ggplot(misclass_rates_logistic, aes(Loan_Risk,Estimated_Misclassification_Rate)) + geom_boxplot(aes(fill = Loan_Risk))

ggplot(correct_default_rates_logistic, aes(Loan_Risk,Estimated_Identitfy_True_Default_Rate)) + geom_boxplot(aes(fill = Loan_Risk))

######################## Discriminant Analysis Hybird Algorithm  ######################
# Segment the data by loan risk category
LC_low = LC_complete %>% filter(loan_risk == "Low")
LC_moderate = LC_complete %>% filter(loan_risk == "Moderate")
LC_high = LC_complete %>% filter(loan_risk == "High")

# Cost ratios by risk level
Cost_ratios_LMH = c(1/60,1/30,1/15)

# Linear Discriminant Analysis with Leave-One-Out-Cross Validattion
# Low risk loans
lda_low = lda(loan_default ~ ., data = LC_low[-c(4,6)], CV = TRUE)

# Add CV posterior estimates to original data
LC_low = cbind(LC_low,
               lda_low$posterior[,1],
               lda_low$posterior[,2])

names(LC_low)[(ncol(LC_low)-1):ncol(LC_low)] = c("PosteriorDensity_0","PosteriorDensity_1")

# Add class proportions 
LC_low$PriorProportion_0 = rep(sum(I(LC_low$loan_default==0))/nrow(LC_low),nrow(LC_low))
LC_low$PriorProportion_1 = 1-LC_low$PriorProportion_0

LC_low = LC_low[!is.na(LC_low$PosteriorDensity_0),]

# Add Prediction based on minimizing cost of misclassification
LC_low$lda_prediction = 
  ifelse((LC_low$PosteriorDensity_1/LC_low$PosteriorDensity_0) >= (Cost_ratios_LMH[1])*(LC_low$PriorProportion_0/LC_low$PriorProportion_1),1,0)

# Moderate risk loans
lda_moderate = lda(loan_default ~ ., data = LC_moderate[-c(4,6)], CV = TRUE)

# Add CV posterior estimates to original data
LC_moderate = cbind(LC_moderate,
               lda_moderate$posterior[,1],
               lda_moderate$posterior[,2])

names(LC_moderate)[(ncol(LC_moderate)-1):ncol(LC_moderate)] = 
  c("PosteriorDensity_0","PosteriorDensity_1")

# Add class proportions 
LC_moderate$PriorProportion_0 = 
  rep(sum(I(LC_moderate$loan_default==0))/nrow(LC_moderate),nrow(LC_moderate))
LC_moderate$PriorProportion_1 = 1-LC_moderate$PriorProportion_0

LC_moderate = LC_moderate[!is.na(LC_moderate$PosteriorDensity_0),]

# Add Prediction based on minimizing cost of misclassification
LC_moderate$lda_prediction = 
  ifelse((LC_moderate$PosteriorDensity_1/LC_moderate$PosteriorDensity_0) >= (Cost_ratios_LMH[2])*(LC_moderate$PriorProportion_0/LC_moderate$PriorProportion_1),1,0)


# High risk loans
lda_high = lda(loan_default ~ ., data = LC_high[-c(4,6)], CV = TRUE)

# Add CV posterior estimates to original data
LC_high = cbind(LC_high,
                    lda_high$posterior[,1],
                    lda_high$posterior[,2])

names(LC_high)[(ncol(LC_high)-1):ncol(LC_high)] = 
  c("PosteriorDensity_0","PosteriorDensity_1")

# Add class proportions 
LC_high$PriorProportion_0 = 
  rep(sum(I(LC_high$loan_default==0))/nrow(LC_high),nrow(LC_high))
LC_high$PriorProportion_1 = 1-LC_high$PriorProportion_0

LC_high = LC_high[!is.na(LC_high$PosteriorDensity_0),]

# Add Prediction based on minimizing cost of misclassification
LC_high$lda_prediction = 
  ifelse((LC_high$PosteriorDensity_1/LC_high$PosteriorDensity_0) >= (Cost_ratios_LMH[3])*(LC_high$PriorProportion_0/LC_high$PriorProportion_1),1,0)


# Combine results
LC_lda = rbind(LC_low,LC_moderate,LC_high)


# Calculate Performance Stats
# Confusion Matrices by Group
All = with(LC_lda,table(loan_default,lda_prediction))
Low = with(LC_low,table(loan_default,lda_prediction))
Moderate = with(LC_moderate,table(loan_default,lda_prediction))
High = with(LC_high,table(loan_default,lda_prediction))


# Misclassification rate by risk
misclass_rates_lda = c(round(100-(100*(sum(diag(All))/sum(All))),1),
                       round(100-(100*(sum(diag(Low))/sum(Low))),1),
                       round(100-(100*(sum(diag(Moderate))/sum(Moderate))),1),
                       round(100-(100*(sum(diag(High))/sum(High))),1))

# True Defaults Correctly Identified
correct_default_rates_lda = c(round((100*(All[2,2]/sum(All[2,]))),1),
                              round((100*(Low[2,2]/sum(Low[2,]))),1),
                              round((100*(Moderate[2,2]/sum(Moderate[2,]))),1),
                              round((100*(High[2,2]/sum(High[2,]))),1))

expected_error_rates_lda = c(round((100*(All[2,1]/sum(All[,1]))),1),
                             round((100*(Low[2,1]/sum(Low[,1]))),1),
                             round((100*(Moderate[2,1]/sum(Moderate[,1]))),1),
                             round((100*(High[2,1]/sum(High[,1]))),1))



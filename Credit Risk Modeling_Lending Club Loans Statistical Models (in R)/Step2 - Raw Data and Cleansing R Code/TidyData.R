library(dplyr, warn.conflicts = FALSE)

# Read in raw data in csv format (adjust path as needed)
LC_data = read.csv("c:/users/dsvancer/desktop/r directory/lending club data/convert to tidy data/LoansRaw.csv",
                   skip = 1, stringsAsFactors = FALSE)

Zip3Rents = read.csv("c:/users/dsvancer/desktop/r directory/lending club data/convert to tidy data/Zip3Rents.csv",
                     skip = 0, stringsAsFactors = FALSE,
                     colClasses = c("character","character","numeric","character",
                                    "numeric","numeric","numeric","numeric",
                                    "numeric"))

# Find and remove leading or trailing whitespace in character vectors
charVars = which(sapply(LC_data,class)=="character")

LC_data[charVars] = lapply(LC_data[charVars],
                          function(col) trimws(col, which = "both"))

# Locate date variables in the raw data set
dates = which(names(LC_data) %in% c("issue_d","earliest_cr_line","last_pymnt_d",
                                    "next_pymnt_d","last_credit_pull_d"))

# Convert date character variables into R date objects
# Raw date variables are in "Dec-2013" format, add "-01" for day to 
# be able to convert to an R date object
for (i in dates) {
LC_data[,i] = as.Date(paste0(LC_data[,i],"-01"),
                          format = "%b-%Y-%d") }


# Filter data for 36-month loans that have completed their term
LC_data = LC_data %>% filter(issue_d <= as.Date("2013-09-01",format="%Y-%m-%d"),
                             term == "36 months",
                             is.na(next_pymnt_d))

# Locate variables that have only missing values (NAs)
AllMissingVars = which(sapply(LC_data,
                        function(col) {sum(is.na(col))}) == nrow(LC_data))

# Remove these columns from the data
LC_data = LC_data[-AllMissingVars]

# Convert int_rate to numeric
# Remove "%" and then convert to numeric variable
LC_data$int_rate = as.numeric(gsub("%","",LC_data$int_rate))

# Convert emp_length to numeric by recoding values
LC_data$emp_length[LC_data$emp_length == "< 1 year"] = 0
LC_data$emp_length[LC_data$emp_length == "1 year"] = 1
LC_data$emp_length[LC_data$emp_length == "2 years"] = 2
LC_data$emp_length[LC_data$emp_length == "3 years"] = 3
LC_data$emp_length[LC_data$emp_length == "4 years"] = 4
LC_data$emp_length[LC_data$emp_length == "5 years"] = 5
LC_data$emp_length[LC_data$emp_length == "6 years"] = 6
LC_data$emp_length[LC_data$emp_length == "7 years"] = 7
LC_data$emp_length[LC_data$emp_length == "8 years"] = 8
LC_data$emp_length[LC_data$emp_length == "9 years"] = 9
LC_data$emp_length[LC_data$emp_length == "10+ years"] = 10
LC_data$emp_length[LC_data$emp_length == "n/a"] = NA

LC_data$emp_length = as.numeric(LC_data$emp_length)

# Recode home_ownership
LC_data$home_ownership[LC_data$home_ownership == "MORTGAGE"] = "Mortgage"
LC_data$home_ownership[LC_data$home_ownership == "RENT"] = "Rent"
LC_data$home_ownership[LC_data$home_ownership == "OWN"] = "Own"
LC_data$home_ownership[LC_data$home_ownership %in% c("OTHER","NONE")] = "Other"

LC_data$home_ownership = factor(LC_data$home_ownership,
                                levels = c("Mortgage","Rent","Own","Other"),
                                labels = c("Mortgage","Rent","Own","Other"))

# Convert loan status to factor
LC_data$loan_status = factor(LC_data$loan_status,
                             levels = c("Fully Paid","Charged Off"),
                             labels = c("Fully Paid", "Charged Off"))

# Remove % from revol_util and convert to numeric
LC_data$revol_util = as.numeric(gsub("%","",LC_data$revol_util))

# Remove other unnecessary variables
# Most of these would not be available at the time
# of loan approval. Our model is aimed at prediciting
# loan default with information available for new loans
LC_data = 
LC_data[-which(names(LC_data) %in% c("term","url","desc","purpose","sub_grade",
        "title","policy_code","application_type","acc_now_delinq","chargeoff_within_12_mths",
        "bc_open_to_buy","collection_recovery_fee","collections_12_mths_ex_med","delinq_amnt",
        "initial_list_status","last_credit_pull_d","last_pymnt_amnt","last_pymnt_d",
        "mo_sin_old_rev_tl_op","mo_sin_old_il_acct","mo_sin_rcnt_rev_tl_op",
        "mo_sin_rcnt_tl","mths_since_last_delinq","mths_since_last_major_derog",
        "mths_since_last_record","mths_since_recent_inq","mths_since_recent_revol_delinq",
        "num_tl_120dpd_2m","num_tl_30dpd","num_tl_90g_dpd_24m","num_tl_op_past_12m",
        "out_prncp","out_prncp_inv","pymnt_plan","recoveries","total_pymnt",
        "total_pymnt_inv","total_rec_int","total_rec_late_fee","total_rec_prncp",
        "funded_amnt","funded_amnt_inv","tot_coll_amt","total_acc","tot_cur_bal",
        "total_rev_hi_lim","acc_open_past_24mths","avg_cur_bal","num_bc_sats","num_bc_tl",
        "tot_hi_cred_lim","total_bal_ex_mort","total_bc_limit","total_il_high_credit_limit",
        "mths_since_recent_bc","mths_since_recent_bc_dlq","verification_status","num_sats"))]



# Derived Variables

# Indicator for loan default status
LC_data$loan_default = as.numeric(I(LC_data$loan_status == "Charged Off"))

# Loan amount as a percent of annual income
LC_data$pct_loan_income = LC_data$loan_amnt/LC_data$annual_inc

# Months since first credit line
LC_data$months_since_first_credit =  
  round(as.numeric(LC_data$issue_d - LC_data$earliest_cr_line)/30.42,0)

# Loan risk (A-B) low, (C-D) moderate, (E-F) high
LC_data$loan_risk[LC_data$grade %in% c("A","B")] = "Low"
LC_data$loan_risk[LC_data$grade %in% c("C","D")] = "Moderate"
LC_data$loan_risk[LC_data$grade %in% c("E","F","G")] = "High"

LC_data$loan_risk = factor(LC_data$loan_risk,
                           levels = c("Low","Moderate","High"),
                           labels = c("Low","Moderate","High"))

# Left Join Zip3 data to LC_data
LC_data = left_join(LC_data,Zip3Rents[c("Zip3","StateCode","CountyID","CountyName","Latitude",
                                        "Longitude","Rent50_2","pop2000")],
                    by = c("zip_code" = "Zip3"))

# Update variable names to match data set naming convention
names(LC_data)[40:46] = c("state_zip_code","county_id","county_name","county_latitude","county_longitude","median_rent_2_bed","population_2000")


# Adjusted Income (Annual income - median rent * 12)
LC_data$adjusted_annual_inc = LC_data$annual_inc - (12*LC_data$median_rent_2_bed)

# Loan amount as a pct of adjusted (for rent) annual income
LC_data$pct_loan_adj_income = LC_data$loan_amnt/LC_data$adjusted_annual_inc


# Arrange Variables
LC_data = LC_data[c("id","member_id","issue_d","loan_status","loan_default",
                    "loan_amnt","int_rate","installment","grade","loan_risk",
                    "emp_title","emp_length","home_ownership","annual_inc",
                    "adjusted_annual_inc","pct_loan_income","pct_loan_adj_income",
                    "zip_code","addr_state","state_zip_code","county_id",
                    "county_name","county_latitude","county_longitude",
                    "median_rent_2_bed","population_2000","dti","delinq_2yrs",
                    "earliest_cr_line","months_since_first_credit","inq_last_6mths",
                    "open_acc","pub_rec","revol_bal","revol_util","bc_util","mort_acc",
                    "num_accts_ever_120_pd","num_actv_bc_tl","num_actv_rev_tl",
                    "num_il_tl","num_op_rev_tl","num_rev_accts","num_rev_tl_bal_gt_0",
                    "pct_tl_nvr_dlq","percent_bc_gt_75","pub_rec_bankruptcies","tax_liens")]

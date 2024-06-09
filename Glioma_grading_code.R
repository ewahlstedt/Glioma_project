# Libraries used
if (!require(caret)) install.packages('caret')
library(caret)

if (!require(kableExtra)) install.packages('kableExtra')
library(kableExtra) # Guidance taken from here: https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_pdf.pdf

if (!require(knitr)) install.packages('knitr')
library(knitr)

if (!require(rpart)) install.packages('rpart')
library(rpart)

if (!require(rpart.plot)) install.packages('rpart.plot')
library(rpart.plot) # Guidance taken from here: http://www.milbo.org/doc/prp.pdf

if (!require(stringr)) install.packages('stringr')
library(stringr)

if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(xgboost)) install.packages('xgboost')
library(xgboost)


# Download the dataset - method used in Movielens project 

# The dataset used for this report (Glioma Grading Clinical and Mutation Features) can be found here:

#Website and more info on dataset: https://archive.ics.uci.edu/dataset/759/glioma+grading+clinical+and+mutation+features+dataset

#Research paper citation: Tasci, E.; Zhuge, Y.; Kaur, H.; Camphausen, K.; Krauze, A.V.
# Hierarchical Voting-Based Feature Selection and Ensemble Learning Model Scheme for Glioma Grading with Clinical and Molecular Characteristics. 
#Int. J. Mol. Sci. 2022, 23, 14155.
#https://www.semanticscholar.org/reader/992bf4c0b92ef251644ac2854dd1baacd7e42dc5

webfiles <- "glioma+grading+clinical+and+mutation+features+dataset.zip" # Name the file for the environment
if(!file.exists(webfiles))
  download.file("https://archive.ics.uci.edu/static/public/759/glioma+grading+clinical+and+mutation+features+dataset.zip", webfiles)

# There are two files - one where each column has actual data recorded (in this case e.g. female or male, and NOT_MUTATED or MUTATED), and one with binary categorical values encoded as 1 or 0.

info_with_grade <- "TCGA_InfoWithGrade.csv" # File where some columns have been removed and binary variables have been encoded 1 or 0
if(!file.exists(info_with_grade))
  unzip(webfiles, info_with_grade)

mutations_all <- "TCGA_GBM_LGG_Mutations_all.csv" # Original dataset
if(!file.exists(mutations_all))
  unzip(webfiles, mutations_all)

minimal_dataset <- read.csv(info_with_grade)
all_info_dataset <- read.csv(mutations_all)

  
# Create table to reference changes that were made by dataset creators when preprocessing data
table_variables_reference <- data.frame(Variable = c(rep("Grade", 2), "Project", "Case_ID", rep("Gender", 2), "Age_at_diagnosis", "Primary_Diagnosis", rep("Race", 4), rep("All columns representing molecular features", 2)),
                                        Key = c("0 = LGG", "1 = GBM", "Removed", "Removed", "0 = Male", "1 = Female", "Converted from string to continuous value", "Removed", "0 = white", "1 = black or african American", "2 = asian", "3 = american indian or alaska native", "0 = NOT_MUTATED", "1 = MUTATED"))

kbl(table_variables_reference, booktabs = TRUE, escape = TRUE, caption = "Changes made to variables of original dataset in minimal dataset") %>% 
  kable_styling(latex_options = "HOLD_position") %>% 
  collapse_rows(columns = 1, latex_hline = "major", row_group_label_position = "first")


#### Data exploration ####

# Create table to show distribution of gender
table_summary_gender <- minimal_dataset %>% 
  mutate(Gender = ifelse(Gender == 0, "Male", "Female")) %>%
  group_by(Gender) %>% 
  summarize(Count = n()) %>% 
  mutate(Percent = round(Count/sum(Count)*100, 2))

kbl(table_summary_gender, booktabs = TRUE, caption = "Distribution of gender in the minimal dataset") %>% 
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")

# Create table to show distribution of race
table_summary_race <- minimal_dataset %>% 
  mutate(Race = str_replace_all(Race, c("0" = "White", "1" = "Black or African American", "2" = "Asian", "3" = "American Indian or Alaska Native"))) %>%
  group_by(Race) %>% 
  summarize(Count = n()) %>% 
  mutate(Percent = round(Count/sum(Count)*100, 2))

kbl(table_summary_race, booktabs = TRUE, caption = "Distribution of race in the minimal dataset") %>% 
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")


# Create table to show distribution of grade
table_summary_grade <- minimal_dataset %>% 
  mutate(Grade = ifelse(Grade == 0, "LGG", "GBM")) %>% 
  group_by(Grade) %>% 
  summarize(Count = n()) %>% 
  mutate(Percent = round(Count/sum(Count)*100, 2))

kbl(table_summary_grade, booktabs = TRUE, caption = "Distribution of grade in the minimal dataset") %>% 
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")


# Create table to show distribution of age at diagnosis
table_summary_age <- minimal_dataset %>% 
  summarize(Mean = round(mean(Age_at_diagnosis), 0), Median = round(median(Age_at_diagnosis), 0), Youngest = round(min(Age_at_diagnosis), 0), Oldest = round(max(Age_at_diagnosis), 0))

table_summary_age_long <- pivot_longer(table_summary_age, cols = everything(), names_to = "Statistic", values_to = "Value")

kbl(table_summary_age_long, booktabs = TRUE, caption = "Distribution of age at diagnosis in the minimal dataset") %>% 
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")


# Plot distribution of age
minimal_dataset %>% ggplot(aes(x = Age_at_diagnosis)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 1, fill = "grey") + 
  geom_density(colour = "violet", linewidth = 1) + 
  theme_bw() +
  labs(x = "Age at diagnosis", y = "Density", title = "Distribution of Age at Diagnosis") + 
  scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 70, 80, 90, 100))

# Create table of molecular features for plotting
table_mol_feats <- minimal_dataset %>% mutate(Grade = ifelse(Grade == 0, "LGG", "GBM")) %>% 
  pivot_longer(cols = 5:ncol(minimal_dataset), names_to = "mol_feats", values_to = "mol_feats_y_n") %>% # Create column with feature names and column indicating if mutated (1) or not (0)
  group_by(mol_feats, Grade) %>% 
  summarize(Count = sum(mol_feats_y_n), .groups = "keep")

# Plot distribution of molecular features
table_plot_mol_feats <- table_mol_feats %>% 
  group_by(mol_feats) %>% 
  summarize(Total = sum(Count))

table_plot_mol_feats %>% 
  ggplot(aes(x = reorder(mol_feats, -Total), y = Total, fill = mol_feats)) +
  geom_col(show.legend = FALSE) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  labs(x = "Molecular Feature", y = "Total Number of Mutation Occurrences", title = "Distribution of Molecular Features")

# Create table for plotting distribution of grade for different genders
table_gender_grade <- minimal_dataset %>% 
  mutate(Gender = ifelse(Gender == 0, "Male", "Female"), Grade = ifelse(Grade == 0, "LGG", "GBM")) %>% 
  group_by(Grade, Gender) %>% 
  summarize(Count = n(), .groups = "keep") %>%
  group_by(Gender) %>%
  summarize(Grade = Grade, Count = Count, Percent = round(Count/sum(Count)*100, 1), .groups = "keep")


# Plot distribution of grade for different genders
table_gender_grade %>% ggplot(aes(x = Grade, y = Count, fill = Grade), show.legend = FALSE) +
  geom_col(show.legend = FALSE) +
  geom_text(data = table_gender_grade, aes(label = paste(Percent, "%"), vjust = -0.5)) + 
  labs(title = "Grade Distribution for Different Genders") +
  theme_bw() +
  facet_grid(. ~ Gender)

# Create table for plotting distribution of grade across different races
table_race_grade <- minimal_dataset %>% 
  mutate(Race = str_replace_all(Race, c("0" = "White", "1" = "Black or African American", "2" = "Asian", "3" = "American Indian or Alaska Native"))) %>%
  mutate(Grade = ifelse(Grade == 0, "LGG", "GBM")) %>%
  group_by(Grade, Race) %>%
  summarize(Count = n(), .groups = "keep") %>%
  group_by(Race) %>%
  summarize(Grade = Grade, Count = Count, Percent = round(Count/sum(Count)*100, 1), .groups = "keep")

# Plot distribution of grade across different races
table_race_grade %>% ggplot(aes(x = Grade, y = Count, fill = Grade), show.legend = FALSE) +
  geom_col(show.legend = FALSE) +
  geom_text(data = table_race_grade, aes(label = paste(Percent, "%"), vjust = -0.5)) + 
  scale_y_log10() + 
  theme_bw() +
  facet_grid(. ~ Race) +
  labs(title = "Grade distribution for different races") 

# Plot differences in mutation between grades for all molecular features
table_mol_feats %>% ggplot(aes(x = Grade, fill = Grade)) +
  geom_col(aes(y = Count), show.legend = FALSE) +
  theme_bw() +
  facet_wrap(~ mol_feats) +
  labs(y = "Number of cases with mutation present", title = "Mutation prescence in different glioma grades by each molecular features") +
  scale_y_log10() # This gives a warning about infinite values and 1 row being removed - this is due to one value being 0 but shows up as such in the plot anyway


# Create table and calculate differences between grades in terms of mutation of molecular features
table_diff_mol_feats <- table_mol_feats %>% 
  pivot_wider(names_from = Grade, values_from = Count) %>% 
  mutate(Diff = abs(LGG - GBM), Pct_diff = round(Diff/(LGG + GBM)*100, 2)) %>% arrange(desc(Pct_diff))

kbl(table_diff_mol_feats, escape = TRUE, booktabs = TRUE, caption = "Differences in presence of mutation of molecular features in cases of LGG and GBM") %>%
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")


# Filter table for mutations that occur at least 100 times
table_diff_mol_feats %>% filter((GBM+LGG) > 100) %>% 
  kbl(booktabs = TRUE, caption = "Differences in presence of mutation of most common molecular features in cases of LGG and GBM") %>% 
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")


# Plot age vs grade
minimal_dataset %>% 
  mutate(Grade = ifelse(Grade == 0, "LGG", "GBM")) %>%
  ggplot(aes(x = Grade, y = Age_at_diagnosis, fill = Grade)) +
  geom_boxplot(show.legend = FALSE) +
  theme_bw() + 
  labs(y = "Age at Diagnosis", title = "Age at diagnosis for LGG and GBM grades") +
  scale_y_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))


#### Method ####

# Check that no CaseID appears twice

table_length_Case_ID <- data.frame(Length_unique = length(unique(all_info_dataset$Case_ID)), # Length of unique Case_IDs
                                   Length_all = length(all_info_dataset$Case_ID)) # Length of whole Case_ID variable

kbl(table_length_Case_ID, booktabs = TRUE, caption = "Comparing length of unique CaseId and all CaseId in unprocessed dataset") %>% kable_styling(latex_options = "HOLD_position")



#### Create train and test sets from the minimal_dataset ####

set.seed(1) # Set seed for reproducibility
minimal_dataset$Grade <- ifelse(minimal_dataset$Grade == 0, "LGG", "GBM") # Change back to abbreviations for interpretability
minimal_dataset$Grade <- factor(minimal_dataset$Grade) # Turn into factor for model building

test_index <- createDataPartition(y = minimal_dataset$Grade, times = 1, p = 0.3, list = FALSE) # Create index for 30 percent of data

train_set <- minimal_dataset[-test_index, ] # Select rows not in test index (70 percent of data)
test_set <- minimal_dataset[test_index, ] # Select rows in test index (30 percent of data)


# Create table for checking distribution of grade variable across train and test sets

train_for_tab <- train_set %>% 
  group_by(Grade) %>% 
  summarize(Count = n()) %>% 
  ungroup() %>% 
  mutate(Set = c(rep("Train", 2)), Percent = round(Count/sum(Count)*100, 1))

test_for_tab <- test_set %>% 
  group_by(Grade) %>% 
  summarize(Count = n()) %>% 
  ungroup() %>% 
  mutate(Set = c(rep("Test", 2)), Percent = round(Count/sum(Count)*100, 1))

full_for_tab <- minimal_dataset %>% 
  group_by(Grade) %>% 
  summarize(Count = n()) %>% 
  ungroup() %>% 
  mutate(Set = c(rep("Full", 2)), Percent = round(Count/sum(Count)*100, 1))

table_train_test <- rbind(full_for_tab, train_for_tab, test_for_tab) %>% 
  relocate(Set, .before = Grade)

kbl(table_train_test, booktabs = TRUE, caption = "Distribution of the Grade variable across training and test sets compared to full dataset") %>% 
  collapse_rows(columns = 1, latex_hline = "major", row_group_label_position = "first") %>%
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")

rm(test_index, train_for_tab, test_for_tab, full_for_tab)


#### Model testing ####

### Decision Tree model with rpart with all features ###

# Set up the model 
set.seed(1)
model_rpart <- train(Grade ~ . , 
                     method = "rpart", 
                     data = train_set,
                     tuneGrid = data.frame(cp = seq(0, 0.1, len = 50)), # Grid search for best cp
                     trControl = trainControl(method = "cv", number = 10) # 10-fold cross validation
)

# Predict on test set
pred_rpart <- predict(model_rpart, test_set)

# Save accuracy and feature importance
acc_rpart <- mean(pred_rpart == test_set$Grade) # Find model accuracy
imp_rpart <- varImp(model_rpart, scale = TRUE) # Find importance of different features

# Plot the decision tree
rpart.plot(model_rpart$finalModel, extra = "auto", box.palette = "BlGnYl")

# Plot accuracy vs complexity parameter (cp)
plot(model_rpart)

# Accuracy and confusion matrix tables for the rpart model

kbl(acc_rpart, booktabs = TRUE, col.names = NULL) %>% kable_styling(latex_options = "HOLD_position", position = "left")

conf_mat_rpart <- confusionMatrix(pred_rpart, test_set$Grade)

kbl(conf_mat_rpart$table, booktabs = F) %>% 
  kable_styling(latex_options = "HOLD_position", position = "left") %>%
  add_header_above(c("Prediction", "Reference" = 2))

# Plot feature importance rpart model
plot(imp_rpart, main = "Importance of Features")


### Random forest model with all predictors ###

# Set up random forest model
set.seed(1)
model_rf <- train(Grade ~ . , 
                  method = "rf", 
                  data = train_set,
                  tuneGrid = data.frame(mtry = seq(1, 5, 1)), # Grid search for best mtry
                  trControl = trainControl(method = "cv", number = 10) # 10-fold cross validation
)

# Predict on test set
pred_rf <- (predict(model_rf, test_set))

# Save accuracy and feature importance
acc_rf <-  mean(pred_rf == test_set$Grade) # Find model accuracy
imp_rf <- varImp(model_rf, scale = TRUE) # Find importance of different features

# Plot random forest model
plot(model_rf)

# Accuracy and confusion matrix tables for the random forest model

kbl(acc_rf, booktabs = TRUE, col.names = NULL) %>% kable_styling(latex_options = "HOLD_position", position = "left")

conf_mat_rf <- confusionMatrix(pred_rf, test_set$Grade)

kbl(conf_mat_rf$table, booktabs = F) %>% 
  kable_styling(latex_options = "HOLD_position", position = "left") %>% 
  add_header_above(c("Prediction", "Reference" = 2))

# Plot feature importance randeom forest model
plot(imp_rf)


### xgboost model with all predictors ###

# Define tuning grid
tune_grid <- expand.grid(nrounds = 200,
                         max_depth = c(3, 6, 9), # Depth of a tree - high values may lead to overfitting
                         eta = c(0.01, 0.1, 0.3), # Controls contribution of each tree to overall ensemble - low values require more rounds to achieve results
                         gamma = c(0, 0.1, 0.2), # Minimum loss required to create new node
                         colsample_bytree = c(0.7, 0.8, 0.9), # Fraction of features to consider when building tree - low value can reduce overfitting
                         min_child_weight = c(1, 2), # Controls if node is split based on size - high values makes model more conservative
                         subsample = c(0.7, 0.8, 0.9) # Fraction of rows to consider when building tree - low value can reduce overfitting
)

# Define train control options
train_control <- trainControl(method = "cv", number = 10, # 10-fold cross-validation
                              verboseIter = FALSE, # Don't print progress to console
                              classProbs = TRUE, # Compute class probabilities
                              summaryFunction = twoClassSummary # As target variable is binary
)
# Set up xgboost model
set.seed(1)
model_xgboost <- train(Grade ~ . , 
                       method = "xgbTree", 
                       data = train_set,
                       tuneGrid = tune_grid, 
                       trControl = train_control
)

# Predict on test set
pred_xgboost <- predict(model_xgboost, test_set)

# Save accuracy and feature importance
acc_xgboost <- mean(pred_xgboost == test_set$Grade) # Find model accuracy
imp_xgboost <- varImp(model_xgboost, scale = TRUE)  # Find importance of different features


# Accuracy and confusion matrix tables for the xgboost model

kbl(acc_xgboost, booktabs = TRUE,col.names = NULL, position = "left")

conf_mat_xgboost <- confusionMatrix(pred_xgboost, test_set$Grade)

kbl(conf_mat_xgboost$table, booktabs = F) %>% 
  kable_styling(latex_options = "HOLD_position", position = "left") %>%
  add_header_above(c("Prediction", "Reference" = 2))


# Plot feature importance xgboost model
plot(imp_xgboost)


#### Feature selection ####

# Check near-zero variance
  
nzv <- nearZeroVar(minimal_dataset)

minimal_dataset[nzv] %>% colnames() %>% kbl(booktabs = TRUE, col.names = NULL, caption = "Features with near-zero variance") %>% kable_styling(latex_options = "HOLD_position")


# Create tables with top five most important features from each model

top_five_rpart <- top_n(data.frame(imp_rpart$importance), 5)
top_five_rpart <- arrange(top_five_rpart, desc(Overall))

top_five_rf <- top_n(data.frame(imp_rf$importance), 5)
top_five_rf <- arrange(top_five_rf, desc(Overall))

top_five_xgboost <- top_n(data.frame(imp_xgboost$importance), 5)
top_five_xgboost <- arrange(top_five_xgboost, desc(Overall))

trpart <- kbl(top_five_rpart, booktabs = TRUE, caption = "rpart model top five features") %>% kable_styling(latex_options = "HOLD_position")

trf <- kbl(top_five_rf, booktabs = TRUE, caption = "random forest model top five features") %>% kable_styling(latex_options = "HOLD_position")

txgboost <- kbl(top_five_xgboost, booktabs = TRUE, caption = "xgboost model top five features") %>% kable_styling(latex_options = "HOLD_position")


#### Re-running models with fewer features ####

### rpart model with selected features ###
# Set up the model 
set.seed(1)
model_rpart_top <- train(Grade ~ IDH1 + Age_at_diagnosis + CIC + ATRX + PTEN + IDH2 + NF1, 
                     method = "rpart", 
                     data = train_set,
                     tuneGrid = data.frame(cp = seq(0, 0.1, len = 50)), # Grid search for best cp
                     trControl = trainControl(method = "cv", number = 10) # 10-fold cross validation
                     )

# Predict on test set
pred_rpart_top <- predict(model_rpart_top, test_set)

# Save accuracy
acc_rpart_top <- mean(pred_rpart_top == test_set$Grade) # Find model accuracy


### Random forest model with selected features ####

# Set up random forest model
set.seed(1)
model_rf_top <- train(Grade ~ IDH1 + Age_at_diagnosis + CIC + ATRX + PTEN + IDH2 + NF1, 
                     method = "rf", 
                     data = train_set,
                     tuneGrid = data.frame(mtry = seq(1, 5, 1)), # Grid search for best mtry
                     trControl = trainControl(method = "cv", number = 10) # 10-fold cross validation
                     )

# Predict on test set
pred_rf_top <- (predict(model_rf_top, test_set))

# Save accuracy
acc_rf_top <- mean(pred_rf_top == test_set$Grade) # Find model accuracy


### xgboost model with selected features ###

# Define tuning grid

tune_grid <- expand.grid(nrounds = 200,
                         max_depth = c(3, 6, 9), # Depth of a tree - high values may lead to overfitting
                         eta = c(0.01, 0.1, 0.3), # Controls contribution of each tree to overall ensemble - low values require more rounds to achieve results
                         gamma = c(0, 0.1, 0.2), # Minimum loss required to create new node
                         colsample_bytree = c(0.7, 0.8, 0.9), # Fraction of features to consider when building tree - low value can reduce overfitting
                         min_child_weight = c(1, 2), # Controls if node is split based on size - high values makes model more conservative
                         subsample = c(0.7, 0.8, 0.9) # Fraction of rows to consider when building tree - low value can reduce overfitting
)

# Define train control options
train_control <- trainControl(method = "cv", number = 10, # 10-fold cross-validation
                              verboseIter = FALSE, # Don't print progress to console
                              classProbs = TRUE, # Compute class probabilities
                              summaryFunction = twoClassSummary # As target variable is binary
)
# Build the model
set.seed(1)
model_xgboost_top <- train(Grade ~ IDH1 + Age_at_diagnosis + CIC + ATRX + PTEN + IDH2 + NF1, 
                           method = "xgbTree", 
                           data = train_set,
                           tuneGrid = tune_grid, 
                           trControl = train_control
)

# Predict on test set
pred_xgboost_top <- predict(model_xgboost_top, test_set)

# Save accuracy
acc_xgboost_top <- mean(pred_xgboost_top == test_set$Grade) # Find model accuracy


# Table of accuracies so far
table_all_accs <- data.frame(Model = c("rpart - all features", "rpart - top features", "random forest - all features", "random forest - top features", "xgboost - all features", "xgboost - top features"), 
                             Accuracy = c(acc_rpart, acc_rpart_top, acc_rf, acc_rf_top, acc_xgboost, acc_xgboost_top))
table_all_accs <- arrange(table_all_accs, desc(Accuracy))

kbl(table_all_accs, booktabs = TRUE, caption = "Accuracy of all models tested so far") %>% 
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")


### Ensemble model ###

# Majority voting ensemble
ensemble_vote <- cbind(rpart = pred_rpart_top == "LGG",  rf = pred_rf_top == "LGG", xgboost = pred_xgboost_top == "LGG")

# Predict on test set
pred_ensemble_vote <- ifelse(rowMeans(ensemble_vote) >= 0.5, "LGG", "GBM")

# Compute and save accuracy
acc_ensemble_vote <- mean(pred_ensemble_vote == test_set$Grade)

  
#### Results ####
  
# Table of accuracies for all models including ensemble
table_all_accs <- data.frame(Model = c("rpart - all features", "rpart - top features", "random forest - all features", "random forest - top features", "xgboost - all features", "xgboost - top features", "Ensemble - voting with top features"), 
                             Accuracy = c(acc_rpart, acc_rpart_top, acc_rf, acc_rf_top, acc_xgboost, acc_xgboost_top, acc_ensemble_vote))
table_all_accs <- arrange(table_all_accs, desc(Accuracy))

kbl(table_all_accs, booktabs = TRUE, caption = "Accuracy of all models tested including ensemble model") %>% 
  kable_styling(full_width = TRUE, latex_options = "HOLD_position") %>% 
  column_spec(1, width = "7cm")

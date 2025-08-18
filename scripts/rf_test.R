#!/usr/bin/env Rscript

# Script header ----
# Programmed by: Kairi Tanaka
# Programmed on: 05-07-2025
# Programmed to: apply random forest classifier for sewage analysis


# 1) Load required packages
if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}
if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(randomForest)
library(caret)
library(ggplot2)

# 2) Parameters: edit these
input_file   <- 
  sample_data(ps_gm) |> 
  data.frame() |> 
  select(Community2, Covid) |> 
  merge(
    otu_table(ps_gm) |> 
      t(),
    by = 'row.names'
  ) |> 
  dplyr::rename(SampleID = 'Row.names') |> 
  filter(Covid !='NA')
# input_file |> View
group_column <- "Community2"
# group_column <- "Covid"
train_frac   <- 0.7
ntree        <- 500

# Read data ----
dat <- input_file

# Clean colnames
colnames(dat) <- make.names(colnames(dat))

# Set SampleID as rownames and drop it
if ("SampleID" %in% colnames(dat)) {
  rownames(dat) <- dat$SampleID
  dat$SampleID  <- NULL
}

# Extract response and features ----
dat[[group_column]] <- as.factor(dat[[group_column]])
features <- setdiff(colnames(dat), group_column)
X <- dat[, features]
y <- dat[[group_column]]

# Split into training and test sets ----
set.seed(123)  # for reproducibility
train_idx <- createDataPartition(y, p = train_frac, list = FALSE)
train_dat <- dat[train_idx, ]
test_dat  <- dat[-train_idx, ]

# Train the random forest model ----
rf_model <- randomForest(
  # formula = as.formula(paste(group_column, "~ .")),
  formula = Community2 ~ Covid + .,
  # formula = Covid ~ Community2 + .,
  data    = train_dat,
  ntree   = ntree,
  importance = TRUE,
  proximity = TRUE
)
print(rf_model)

# Quick tuning ----
best <- tuneRF(
  x = train_dat[, features],
  y = train_dat$Community2,
  ntreeTry = 500,
  stepFactor = 1.5,
  improve = 0.01,
  trace = TRUE,
  plot = TRUE
)

# CV grid search ----
library(caret)
ctrl <- trainControl(
  method      = "repeatedcv",
  number      = 5,
  repeats     = 3,
  # classProbs  = FALSE
  classProbs  = TRUE, # if using AUC
  summaryFunction = multiClassSummary
)

tunegrid <- expand.grid(
  mtry = seq(5, floor(sqrt(ncol(train_dat)-1)*3), by=2)
)

rf_caret <- train(
  Community2 ~ Covid + .,
  data      = train_dat,
  method    = "rf",
  metric    = "Accuracy",
  tuneGrid  = tunegrid,
  trControl = ctrl,
  ntree     = 500
)

plot(rf_caret)  # ROC or accuracy vs mtry

# Evaluate on test set ----
pred <- predict(rf_model, newdata = test_dat)
conf <- confusionMatrix(pred, test_dat[[group_column]])
cat("\nConfusion Matrix:\n")
print(conf$table)
cat("\nOverall Statistics:\n")
print(conf$overall)

# OOB error
plot(rf_model)

# Partial dependence plot
partialPlot(rf_model, pred.data=train_dat, x.var="Thauera", which.class="Minority")

# MDS on proximity matrix
MDSplot(rf_model, train_dat$Community2, palette=RColorBrewer::brewer.pal(3,'Dark2'), pch=as.numeric(train_dat$Community2))

# PCA on proximity matrix
pca <- princomp(rf_model$proximity)
factoextra::fviz_eig(pca)
pca$scores |> 
  merge(sample_data(ps_filt), by = 'row.names') |> 
  ggplot(mapping=aes(x=Comp.1, y=Comp.2)) +
  geom_point(mapping=aes(shape=Community2, color=Community2))

# Plot variable importance ----
# Top 20 taxa
imp <- importance(rf_model, type = 1)
imp_df <- data.frame(
  Taxa       = rownames(imp),
  MeanDecreaseAccuracy = imp[, 1],
  row.names  = NULL
)
imp_df <- imp_df[order(imp_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
top20  <- head(imp_df, 20)

ggplot(top20, aes(x = reorder(Taxa, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
# ggplot(imp_df, aes(x = reorder(Taxa, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 20 Important Taxa (Mean Decrease in Accuracy)",
    x     = "",
    y     = "Mean Decrease Accuracy"
  ) +
  theme_cowplot()

# 9) Save the model and importance table
# saveRDS(rf_model, file = "rf_taxa_model.rds")
# write.csv(imp_df, file = "taxa_importance.csv", row.names = FALSE)


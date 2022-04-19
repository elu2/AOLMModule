# Path of RDS file which xy data is stored in
xy_rds <- ""
# Path of directory to write alpha grids and test indices to
grid_path <- "./grids/"
if (!(file.exists(grid_path))){dir.create(grid_path)}
# Path of csv to write validation predictions to
val_csv <- "./validationOut.csv"

# For frequency validation
n_iter <- 100
freq_val_csv <- "./freqValOut.csv"

# glmnet "family" argument
family <- "binomial"
# glmnet "type.measure" argument
type.measure="class"

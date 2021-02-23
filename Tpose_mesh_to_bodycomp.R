# This script allows you to use the pre-existing PCA model to generate PC weights using new meshes and generate 
# body composition values

# load dependencies
library(Rvcg)

# Set working directory - make sure meshes are in here
setwd('C:/Users/...')

# load xyz coordinates from PLY file to vector
ply_to_vec <- function(ply_file) {
  cat("Loading:", ply_file, "\n")
  pts = vcgPlyRead(ply_file)
  xyz = pts$vb[1:3,]
  vec = as.vector(xyz)
  return(vec)
}

# load all PLY files in directory
files = list.files(path = getwd(), pattern = ".ply")

# parse by sex if your meshes are not already separated, go to the sex_id.csv to input the Fit3D package ID and sex 
df_sex_id = read.csv('sex_id.csv')

files_m = c()
files_f = c()
ids_m = c()
ids_f = c()

for (file in files) {
  if (startsWith(file, "pkg")) {
    # file is named with package ID
    package_id = substr(file, 0, 36)
    rows = df_sex_id[df_sex_id$Scan_Package_ID == package_id,]
  } else if (grepl("[0-9]{2}ADL[0-9]{4}", file)) {
    # file is named with SubjectID
    subject_id = substr(file, 0, 9)
    rows = df_sex_id[df_sex_id$SubjectID == subject_id,]
  } else {
    # file has name of unknown format
    cat("Unknown filename: ", file, "\n")
  }
  row = rows[1,]
  sex = row$SEX
  su_id = as.character(row$SubjectID)
  if (length(sex) < 1) {
    cat("Could not find match for file: ", file, "\n")
  } else {
    if (is.na(sex)) {
      cat("Missing sex: ", file, "\n")
    } else if (sex == "M") {
      files_m = c(files_m, file)
      ids_m = c(ids_m, su_id)
    } else if (sex == "F") {
      files_f = c(files_f, file)
      ids_f = c(ids_f, su_id)
    } else {
      cat("Unknown sex: ", sex, " for file ", file, "\n")
    }
  }
}

# upload PCA model
load('pca_model_centered_Tpose_training_121919.Rdata')

# read ply
male_ply = do.call(rbind, lapply(files_m, ply_to_vec))
female_ply = do.call(rbind, lapply(files_f, ply_to_vec))

# use existing pca model to get pc weights
m_pca = scale(male_ply, pca_m$center, pca_m$scale) %*% pca_m$rotation
f_pca = scale(female_ply, pca_f$center, pca_f$scale) %*% pca_f$rotation

# save weight matrices
pcw_m = cbind(meshfile_name = files_m, m_pca)
pcw_f = cbind(meshfile_name = files_f, f_pca)

# turn matrices into datafram
pcw_m = as.data.frame(pcw_m[,1:16])
pcw_f = as.data.frame(pcw_f[,1:16])

# T-pose body composition
pcw_f$Tpose_fatmass = 24.99 + .11*pcw_f$PC1 + 2.04*pcw_f$PC2 - .92*pcw_f$PC3 -.56*pcw_f$PC4 +.53*pcw_f$PC5 - 1.24*pcw_f$PC7 + 
  1.15*pcw_f$PC9 + 1.93*pcw_f$PC10 + 1.93*pcw_f$PC12 - 2.49*pcw_f$PC13 + 1.53*pcw_f$PC15
pcw_m$Tpose_fatmass = 19.66-(.15*pcw_m$PC1)+(1.59*pcw_m$PC2)-(1.56*pcw_m$PC3)+(1.12*pcw_m$PC4)+(.8*pcw_m$PC5)+(2.14*pcw_m$PC7)-(1.58*pcw_m$PC8)+(1.75*pcw_m$PC10)+(1.32*pcw_m$PC13)

# save
write.csv(pcw_m, "male.csv", row.names = FALSE)
write.csv(pcw_f, "female.csv", row.names = FALSE)

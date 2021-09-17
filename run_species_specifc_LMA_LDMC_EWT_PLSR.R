####################################################################################################
#
#
#
#   This R script is provided as supplementary material of Helsen & Bassi et al. 2021 Ecological Indicators.
#
#   The script illustrate how to estimate LMA, LDMC and EWT using the publication dataset 
#   and the provided species specific model for Rosa rugosa, Rubus caesius, Jacobea vulgaris 
#   and Hieracium umbellatum. 
#   
#
#
#    Notes:
#     * The script are re-adapted from the one published by Serbin et al. 2019 New Phytologist
#     * The author notes the code is not the most elegant or clean, but is functional 
#     * Questions, comments, or concerns can be sent to leonardo.bassii@gmail.com
#     * Code is provided under GNU General Public License v3.0 
#
#
#    --- Last updated:  30.08.2021 By Leonardo Bassi <leonardo.bassii@gmail.com>
####################################################################################################

#---------------- Close all devices, delete all variables and prepare all libraries ----------------------

rm(list=ls())   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

### get all required libraries to run script
list.of.packages.CRAN <- c("ggplot2","remotes","devtools","readr","RCurl","httr","pls","dplyr",
                           "reshape2","here","plotrix","gridExtra","scales","knitr", "tidyr","cowplot")  

### check for dependencies and install if needed
new.packages.CRAN <- list.of.packages.CRAN[!(list.of.packages.CRAN %in% installed.packages()[,"Package"])]
if(length(new.packages.CRAN)) install.packages(new.packages.CRAN)

# Check for github 'spectratrait' library and install if needed 
# Note: * you should pay careful attention at this stage to any 
#       *R messages in your terminal alerting 
#       * you that you need to update existing or install new R packages

if(!"spectratrait" %in% installed.packages()[,"Package"]) 
  devtools::install_github(repo = "TESTgroup-BNL/PLSR_for_plant_trait_prediction", dependencies=TRUE) 

### Load libraries 
library(dplyr)
library(ggplot2)


rm(list.of.packages.CRAN,new.packages.CRAN)
#' 
#' 
#-----------------Set working directory (scratch space)-----------------------------------------------

wd <- 'scratch'  #CHANGE THIS TO YOUR PREFERRED DIRECTORY 
require("knitr")
if (! file.exists(wd)) dir.create(file.path("~",wd),recursive=TRUE, showWarnings = FALSE)
setwd(file.path("~",wd)) # set working directory
getwd()  # check wd

#' 
#' 
#-----------------Grab PLSR Coefficients and plot list from GitHub and data from EcoSIS-------------------------------

### Grab Coefficient 

git_repo <- "https://raw.githubusercontent.com/LeonardoUU/Helsen-Bassi_et_al_2021_Ecological_Indicators/main/"
githubURL_coeffs <- paste0(git_repo,"Coefficents_Species-specifc_RMSEP_PLSR_model_LMA_LDMC_EWT.csv")
PLSR.coeffs <- spectratrait::source_GitHubData(githubURL_coeffs)


### Grab list of plot used for validation (exclude data point used to build models)
githubURL_val_plot_list <- paste0(git_repo,"Validation_data_ids.csv")
validation_plot_list<- spectratrait::source_GitHubData(githubURL_val_plot_list)

  #read.csv("Validation_data_ids.csv")

rm(githubURL_coeffs,githubURL_val_plot_list,git_repo )

### Grab Example data 

#### URL:  https://ecosis.org/package/e88b832d-d7da-48b5-af59-25a8079a0ab6  (publication Species specifc dataset)

ecosis_id  <- "e88b832d-d7da-48b5-af59-25a8079a0ab6"  #
dat_raw <- spectratrait::get_ecosis_data(ecosis_id = ecosis_id)
trait_spectra_df<- dat_raw[dat_raw$ids %in% as.vector(validation_plot_list$ids), ] # remouve obs. used to build models


rm(ecosis_id, dat_raw, validation_plot_list)

#' 
#' 
# ----------------Plot data -----------------------------------------------------------------

### Calculate means and SD for plotting
waves <- seq(400,2450,1)
colnames(trait_spectra_df[,1:10])

Plot_data<- trait_spectra_df%>%
  dplyr::select( `Latin Genus`, as.character(waves) )%>%  # grab column 
  tidyr::pivot_longer(!`Latin Genus`, names_to = "wavelength", values_to = "Reflectance")%>% # make long format 
  dplyr::group_by(`Latin Genus`, wavelength)%>%  
  dplyr::summarise(mean_ref = mean(Reflectance), 
                   sd_up = mean(Reflectance) + sd(Reflectance),
                   sd_low = mean(Reflectance) - sd(Reflectance))%>% # calculate mean and standard deviation per species per band
  tidyr::pivot_longer(c("mean_ref","sd_up", "sd_low")  , names_to = "Sum_stat", values_to = "Reflectance")%>% # make long format
  dplyr::mutate(wavelength= as.numeric(wavelength), # get wavelength as numeric for plotting
                line_type_legend= substr(Sum_stat, 1,2), ## make variable for legened 
                `Latin Genus`= factor(`Latin Genus`, 
                                        levels = c("Rosa",  
                                                   "Rubus" ,
                                                   "Jacobaea",
                                                   "Hieracium"))) # Reorder species for plotting
  

### make plot 
plot_spectra_summary<- ggplot2::ggplot(data = Plot_data, aes(x= wavelength,y= Reflectance, group= Sum_stat,
                                             linetype=line_type_legend ))+
  geom_line()+
  scale_linetype_manual(name= "Summary spectra",
                        labels=c("Mean","+ / - SD",""),
                        values=c("solid", "dashed"))+
  labs(x = "Wavelength (nm)",
     y = "Reflectance",
     linetype= "Legend") +
  theme(legend.position="top")+
  facet_wrap(~`Latin Genus`, nrow = 2, ncol = 2, scales = "fixed")+
  theme_bw()+
  theme(legend.position="top",text = element_text(size=15))


### save plot 
png(file=file.path("~",wd,'Species_spectra_summary_plot.png'),height=3000,
    width=3900, res=340)

plot_spectra_summary

dev.off()


 rm(plot_spectra_summary,Plot_data)
 
#' 
#' 
## ---------------Applying PLSR model to estimate LMA, LDMC, EWT from spectral observations-------------------------

### setup model
models_intercept<- PLSR.coeffs[1,]
models_coeffs_matrix<- as.matrix(PLSR.coeffs[-1,-1])

rm(PLSR.coeffs)
### Set up data 
spectra_matrix<-as.matrix(trait_spectra_df[,which(names(trait_spectra_df) %in% seq(400,2450,1))],)
row.names(spectra_matrix) <- trait_spectra_df$ids

### Run models to estimate LMA, LDMC and EWT, with each of the species specif model

temp<- spectra_matrix%*% models_coeffs_matrix # multiply spectra and coefficient matrix
model_names<- colnames(temp) # grab models names to add intercept using sapply

# define fun. to add intercpet with sapply
add_intercept_fun<- function(model,intercpet, model_name){
  model[,model_name]+intercpet[,model_name]} 

# add intercept with sapply and 
trait_prediction<- sapply(model_names, function(z) add_intercept_fun(temp,models_intercept,z),
                          USE.NAMES=TRUE)%>%
  as.data.frame()%>%
  tibble::rownames_to_column(., var = "ids")

colnames(trait_prediction)<-sub("Coeffs", "Predicted",colnames(trait_prediction)) # rename columns

### Combine measured traits with estimated triats ###
`%!in%` = Negate(`%in%`)
trait_measured_cross_prediction<-  merge(trait_spectra_df[,which(names(trait_spectra_df) %!in% seq(350,2500,1))],
                                   trait_prediction, by="ids")

rm(models_intercept,models_coeffs_matrix,spectra_matrix,temp)


### Extract trait prediction of each species using model trained with the same species

# define function to grab correct model for correct species and rename variable
extr_dat_fun<- function(data, species){
  
  col_grep_expr<-paste(species, sep = "|", "Latin|ids|cm2|mg")
  self_PLS_predicion<-data[grep(species, data$`Latin Genus`), 
                            grep(col_grep_expr, colnames(data))]
  # rename columns
  colnames(self_PLS_predicion)[grep("^Predicted_LMA",colnames(self_PLS_predicion))] <- 'Predicted_LMA_own_PLSR_model'
  colnames(self_PLS_predicion)[grep("^Predicted_LDMC",colnames(self_PLS_predicion))] <- 'Predicted_LDMC_own_PLSR_model'
  colnames(self_PLS_predicion)[grep("^Predicted_EWT",colnames(self_PLS_predicion))] <- 'Predicted_EWT_own_PLSR_model'
  return(self_PLS_predicion)
  }

# use function to create new dataframe with observed traits and traits predicted with species own model 
species_names_list<- c("Rosa","Rubus","Jacobaea","Hieracium")
names(species_names_list)<-species_names_list

trait_measured_self_prediction<- lapply(species_names_list, function(x)extr_dat_fun(trait_measured_cross_prediction,x) )%>%
   do.call("rbind", .)


### Export dataframe
write.csv(x = trait_measured_cross_prediction, file = file.path("~",wd,"Species_specifc_PLSR_corss_predicted_LMA_LDMC_EWT_data.csv"),
                     row.names = FALSE)

write.csv(x = trait_measured_self_prediction, file = file.path("~",wd,"Species_specifc_PLSR_self_predicted_LMA_LDMC_EWT_data.csv"),
          row.names = FALSE)


#----------------- Get R2 and RMSE  self prediction --------------------------------------------------------

model_fit<- function(data, 
                     species){
  
  sp_data<- data[grep(species, data$`Latin Genus`),]

  rmse_lma<- sqrt(mean(( sp_data$`Leaf mass per area (mg/cm2)` - sp_data$Predicted_LMA_own_PLSR_model)^2))
  rmse_ldmc<- sqrt(mean(( sp_data$`Leaf dry matter content (mg/mg)` - sp_data$Predicted_LDMC_own_PLSR_model)^2))
  rmse_ewt<- sqrt(mean(( sp_data$`EWT (mg/cm2)` - sp_data$Predicted_EWT_own_PLSR_model)^2))
  
  r2_lma<- summary(lm(`Leaf mass per area (mg/cm2)`~Predicted_LMA_own_PLSR_model ,data=sp_data))$r.squared
  r2_ldmc<- summary(lm(`Leaf dry matter content (mg/mg)`~Predicted_LDMC_own_PLSR_model ,data=sp_data))$r.squared
  r2_ewt<- summary(lm(`EWT (mg/cm2)`~Predicted_EWT_own_PLSR_model ,data=sp_data))$r.squared
  
  
  data_output<- data.frame( trait= c("Leaf mass per area (mg/cm2)", "Leaf dry matter content (mg/mg)", "EWT (mg/cm2)"),
                            R2= c(r2_lma,r2_ldmc,r2_ewt ),
                            RMSE = c(rmse_lma,rmse_ldmc,rmse_ewt)
  )
    
   
  }

model_info_df<- lapply(species_names_list, function(x) model_fit(data=trait_measured_self_prediction,
                                                     species=x)) %>%
  dplyr::bind_rows(., .id = "Species")


# export model info 
write.csv(x= model_info_df,file = file.path("~",wd,"sefl_prediction_R2_and_RMSE.csv"))

## ---------------Plot self prediction results-----------------------------------------------------

### Plot LMA 

LMA_species_self_prediction_plot<- ggplot(trait_measured_self_prediction, aes(x= Predicted_LMA_own_PLSR_model, 
                                                          y= `Leaf mass per area (mg/cm2)`,  
                                                           color=`Latin Genus`,shape = `Latin Genus`))+
  geom_point(size=0.5)+
  geom_smooth(method="glm", formula = y ~ x)+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  xlab("Predicted LMA via PLSR")+ 
  ylab("Observed LMA (g/m2)")+
  scale_color_manual(values = c("gold2", "darkred", "dodgerblue4", "darkgreen" ))+
  scale_shape_manual(values=c(16,17,18,15))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() , 
        text=element_text(size=12),
        axis.text.x= element_text(colour="black"),
        axis.text.y= element_text(colour="black")) +
  theme(aspect.ratio=1)+
  scale_y_continuous(expand = c(0, 0), # set the orging of the plot to 0, 0
                     limits = c(0, 15), # set range to be used on the axis lables
                     breaks = seq(0, 15, by = 2))+   # set brakes for axis lables
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15), breaks = seq(0, 15, by = 2))



### plot LMDC 
LDMC_species_self_prediction_plot<- ggplot(trait_measured_self_prediction, aes(x= Predicted_LDMC_own_PLSR_model, 
                                                                   y=  `Leaf dry matter content (mg/mg)`,  
                                                                    color=`Latin Genus`,shape = `Latin Genus`))+
  geom_point(size=0.5)+
  scale_color_manual(values = c("gold2", "darkred", "dodgerblue4", "darkgreen" ))+
  scale_shape_manual(values=c(16,17,18,15))+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  geom_smooth(method="glm", formula = y ~ x)+
  ylab("Observed LDMC (mg/g)")+
  xlab("Predicted LDMC via PLSR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() , 
        text=element_text(size=12),
        axis.text.x= element_text(colour="black"),
        axis.text.y= element_text(colour="black")) +
  theme(aspect.ratio=1)+
  scale_y_continuous(expand = c(0, 0), # set the origin of the plot to 0, 0
                     limits = c(0, 0.58), # set range to be used on the axis lables
                     breaks = seq(0, 0.5, by = 0.1))+   # set brakes for axis lables
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.58), breaks = seq(0, 0.6, by = 0.1))

### plot EWT 

EWT_species_self_prediction_plot<- ggplot(trait_measured_self_prediction, aes(x=Predicted_EWT_own_PLSR_model,
                                                                 y=`EWT (mg/cm2)`, 
                                                                 color=`Latin Genus`,shape = `Latin Genus`))+
  geom_point(size=0.5)+
  scale_color_manual(values = c("gold2", "darkred", "dodgerblue4", "darkgreen" ))+
  scale_shape_manual(values=c(16,17,18,15))+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  geom_smooth(method="glm", formula = y ~ x)+
  ylab("Observed EWT (mg/cm2)")+
  xlab("Predicted EWT via PLSR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() , 
        text=element_text(size=12),
        axis.text.x= element_text(colour="black"),
        axis.text.y= element_text(colour="black")) +
  theme(aspect.ratio=1)+
  scale_y_continuous(expand = c(0, 0), # set the orging of the plot to 0, 0
                     limits = c(0, 43), # set range to be used on the axis lables
                     breaks = seq(0, 40, by = 10))+   # set brakes for axis lables
  scale_x_continuous(expand = c(0, 0), limits = c(0, 43), breaks = seq(0, 40, by = 10))



### combine and save plot

legend <-  cowplot::get_legend(LMA_species_self_prediction_plot+
                       theme(legend.position = "top",
                             legend.spacing.x = unit(0.2, 'cm')))

plot_all<- cowplot::plot_grid(
  LMA_species_self_prediction_plot+theme(legend.position = "none"), 
                              LDMC_species_self_prediction_plot+theme(legend.position = "none"),
                              EWT_species_self_prediction_plot+theme(legend.position = "none"),
                              nrow = 1)

plot_all_legend<- cowplot::plot_grid(legend, plot_all,  nrow=2,
                                     rel_heights = c(0.1, 1 ),
                                     rel_widths = c( 0.8, 1))


png(file = file.path("~",wd,"LMA_LMDC_EWT_validation_plot.png"),
    height=1300,
    width=3900, res=340)
plot_all_legend
dev.off()

rm(list=ls())   # clear workspace



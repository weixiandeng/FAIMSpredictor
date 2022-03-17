
#' @importFrom magrittr "%>%"
##########Define function to calculate peptide properties########

Input_matrix_calculation<-function(Seq,Chg){
  x=data.frame(Sequence=Seq,Charge=Chg)
  Pep_aIndex<-Peptides::aIndex(x$Sequence)
  Pep_predicted_charge<-Peptides::charge(x$Sequence, pH = 2.5, pKscale = "EMBOSS")
  PeP_pI<-Peptides::pI(x$Sequence,pKscale = "EMBOSS")
  Pep_InstabilityIndex<-Peptides::instaIndex(x$Sequence)
  Pep_hydrophobicity<-Peptides::hydrophobicity(x$Sequence,scale = "Eisenberg")
  Pep_comp<-Peptides::aaComp(x$Sequence)
  Pep_comp1 = data.frame(matrix(nrow=nrow(x), ncol=18))
  for (i in 1:nrow(Pep_comp1)){Pep_comp1[i,] = as.vector(Pep_comp[[i]])}
  Pep_MW<-Peptides::mw(x$Sequence)
  Pep_length<-Peptides::lengthpep(x$Sequence)
  #new featrues
  Pep_crucianiPeroperties<-Peptides::crucianiProperties(x$Sequence)
  Pep_crucianiPeroperties1 = data.frame(matrix(nrow=nrow(x),ncol=3))
  for (i in 1:nrow(Pep_crucianiPeroperties1))
  {Pep_crucianiPeroperties1[i,] = as.vector(Pep_crucianiPeroperties[[i]])}
  colnames(Pep_crucianiPeroperties1)<-c("PP1","PP2","PP3")

  Pep_fasgaiVectors<-Peptides::fasgaiVectors(x$Sequence)
  Pep_fasgaiVectors1 = data.frame(matrix(nrow=nrow(x),ncol=6))
  for (i in 1:nrow(Pep_fasgaiVectors1))
  {Pep_fasgaiVectors1[i,] = as.vector(Pep_fasgaiVectors[[i]])}
  colnames(Pep_fasgaiVectors1)<-c("F1","F2","F3","F4","F5","F6")

  Pep_hmoment100<-Peptides::hmoment(x$Sequence,angle=100,window=11)
  Pep_hmoment160<-Peptides::hmoment(x$Sequence,angle=160,window=11)

  Pep_kideraFactors<-Peptides::kideraFactors(x$Sequence)
  Pep_kideraFactors1 = data.frame(matrix(nrow=nrow(x),ncol=10))
  for (i in 1:nrow(Pep_kideraFactors1))
  {Pep_kideraFactors1[i,] = as.vector(Pep_kideraFactors[[i]])}
  colnames(Pep_kideraFactors1)<-c("KF1","KF2","KF3","KF4","KF5","KF6",
                                  "KF7","KF8","KF9","KF10")

  Pep_mswhimScores<-Peptides::mswhimScores(x$Sequence)
  Pep_mswhimScores1 = data.frame(matrix(nrow=nrow(x),ncol=3))
  for (i in 1:nrow(Pep_mswhimScores1))
  {Pep_mswhimScores1[i,] = as.vector(Pep_mswhimScores[[i]])}
  colnames(Pep_mswhimScores1)<-c("MSWHIM1","MSWHIM2","MSWHIM3")

  Pep_stscales<-Peptides::stScales(x$Sequence)
  Pep_stscales1 = data.frame(matrix(nrow=nrow(x),ncol=8))
  for (i in 1:nrow(Pep_stscales1))
  {Pep_stscales1[i,] = as.vector(Pep_stscales[[i]])}
  colnames(Pep_stscales1)<-c("ST1","ST2","ST3","ST4","ST5","ST6",
                             "ST7","ST8")

  Pep_tscales<-Peptides::tScales(x$Sequence)
  Pep_tscales1 = data.frame(matrix(nrow=nrow(x),ncol=5))
  for (i in 1:nrow(Pep_tscales1))
  {Pep_tscales1[i,] = as.vector(Pep_tscales[[i]])}
  colnames(Pep_tscales1)<-c("T1","T2","T3","T4","T5")

  Pep_VHSEScales<-Peptides::vhseScales(x$Sequence)
  Pep_VHSEScales1 = data.frame(matrix(nrow=nrow(x),ncol=8))
  for (i in 1:nrow(Pep_VHSEScales1))
  {Pep_VHSEScales1[i,] = as.vector(Pep_VHSEScales[[i]])}
  colnames(Pep_VHSEScales1)<-c("VHSE1","VHSE2","VHSE3","VHSE4","VHSE5","VHSE6",
                               "VHSE7","VHSE8")

  Pep_zscales<-Peptides::zScales(x$Sequence)
  Pep_zscales1 = data.frame(matrix(nrow=nrow(x),ncol=5))
  for (i in 1:nrow(Pep_zscales1))
  {Pep_zscales1[i,] = as.vector(Pep_zscales[[i]])}
  colnames(Pep_zscales1)<-c("Z1","Z2","Z3","Z4","Z5")

  x_aggregated_physichemistry<- cbind(x,
                                      Pep_aIndex,
                                      Pep_predicted_charge,
                                      PeP_pI,
                                      Pep_InstabilityIndex,
                                      Pep_hydrophobicity,
                                      Pep_MW,
                                      Pep_length,
                                      Pep_comp1,
                                      Pep_crucianiPeroperties1,
                                      Pep_fasgaiVectors1,
                                      Pep_hmoment100,
                                      Pep_hmoment160,
                                      Pep_kideraFactors1,
                                      Pep_mswhimScores1,
                                      Pep_stscales1,
                                      Pep_tscales1,
                                      Pep_VHSEScales1,
                                      Pep_zscales1)
  x_input<-x_aggregated_physichemistry%>%
    dplyr::mutate(Seq_Charge=paste(Sequence,"_",Charge))%>%dplyr::select(Seq_Charge,Charge:Z5)

  return(x_input)
}


#' Predict CV value for a given peptide sequence and Charge state pair.
#'
#'
#' @param Sequence A string contains the sequence you want to predict CV value
#' @param Charge Charge state for the peptide sequence you are trying to predict the CV value
#' @return recomended CV value without empirical information
#' @export
CV_calc_no_empirical=function(Sequence,Charge){
  input=Input_matrix_calculation(Sequence,Charge)
  localH2O <- h2o::h2o.init(ip = 'localhost', port = 54321)
  model_path=getwd()
  imported_model <- h2o::h2o.import_mojo(paste0(model_path,"/model1_no_empirical.zip"))
  write.csv(input,file = paste0("input",".csv"))
  datapath <- paste0("input",".csv")
  data.hex <- h2o::h2o.uploadFile(path = datapath)
  data_viriable<-data.hex[,3:ncol(data.hex)]
  predict<-as.data.frame(h2o::h2o.predict(object=imported_model,data_viriable))
  return(predict)
}



#' Predict CV value for a given peptide sequence and Charge state pair, set a empirical value to provide a more precise prediction.
#'
#'
#' @param Sequence A string contains the sequence you want to predict CV value
#' @param Charge Charge state for the peptide sequence you are trying to predict the CV value
#' @param emprirical_cv CV value obtained from prior experiments where the highest intensity was detected
#' @return recomended CV value with empirical information
#' @export
CV_calc_with_empirical=function(Sequence,Charge,empirical_CV){
  input_no_empirical=Input_matrix_calculation(Sequence,Charge)
  input=input_no_empirical%>%mutate(empirical_value=empirical_CV)
  localH2O <- h2o::h2o.init(ip = 'localhost', port = 54321)
  model_path=getwd()
  imported_model <- h2o::h2o.import_mojo(paste0(model_path,"/model2_empirical.zip"))
  write.csv(input,file = paste0("input",".csv"))
  datapath <- paste0("input",".csv")
  data.hex <- h2o::h2o.uploadFile(path = datapath)
  data_viriable<-data.hex[,3:ncol(data.hex)]
  predict<-as.data.frame(h2o::h2o.predict(object=imported_model,data_viriable))
  return(predict)
}



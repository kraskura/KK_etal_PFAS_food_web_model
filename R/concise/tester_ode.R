
# Liang et al 2022 model: https://doi.org/10.1016/j.scitotenv.2022.156397

# write a function that return uptake rate and depends on whatever variables 
# required to make that calculation (e.g, temperature, salinity, species) 
# then instead of multiplying the constant uptake rate by the water PFAS 
# concentration to get your PFAS_in variable,
# use the return of your function call instead

# *******************************************
# PFOS
# ******************************************
# how much time model over 
t_span <- 100 # days
t_step <- 0.5 # day; every ~ 15 min timestep
seq_step<-seq(1, t_span, by = t_step)

# all assumptions
c_f<-100 # fist starting point?  ng/kg; LOW STARTING POINT
c_w<-24.8 # assumed in water ng/L ; Liang et al 2022 Table S1
c_s<-0.77 # assumed in sediment ng/g
c_d<-300 # ng/kg

# example from Liang et al 2022 PFOS
# uptake rates
k1<-1.595788e+01 # gill; L kg-1 d-1
# Sun et al: k1 depends on membrane-water partition coefficient
# Liang et al: k1 depends on k_1=(E_W × G_V)/W ; gill uptake efficiency E_W=〖(1.85+(155/K_ow)〗^(-1)
# G_V=(980×W^0.65)/D_ox 
kd<-0.05750342 # diet kg kg-1 d-1

# elimination rates; d-1
ke<-0.0042992313 # egestion/feces/waste
k2<-1.517242e-02 # gill
kg<-0.001778949 # growth rate constant
kr<-0.0040006888 # renal
# km<-0.0040006888 # metabolic transformation

# kout<-0.037 # Liang et al 2022 PFOS table 1; total elimination out 

# initialize the outputs lists
results_c_f<-c()
results_c_w<-c()
results_c_s<-c()
results_gill_up<-c()
results_diet_up<-c()
results_total_out<-c()

results_c_f<-append(results_c_f, c_f)
results_c_w<-append(results_c_w, c_w)
results_c_s<-append(results_c_s, c_s)
results_gill_up<-append(results_gill_up, 0)
results_diet_up<-append(results_diet_up, 0)
results_total_out<-append(results_total_out, 0)
    
for(i in 1:length(seq_step)){
  
  # Defining the arrows 
  # Uptake via gills, both in kgPFAS/kgFISH/d
  gill_up<-k1*results_c_w[i]
  diet_up<-kd*c_d
  
  total_up<-(gill_up + diet_up)
  
  # Elimination rates, all in d-1  
  renal_out<-kr*results_c_f[i]
  egestion_out<-kd*results_c_f[i]
  gill_out<-k2*results_c_f[i]
  growth_out<-kg*results_c_f[i]
  
  total_out<-(renal_out + egestion_out + gill_out + growth_out)
  
  # change in fish in respect to time 
  # instantaneous rate of change in PFAS in fish
  # (dC_F)/dt=k_1 C_W+k_D  (P_i C_(D,i))-(k_2+k_E+k_M+k_G)C_F # Liang et al 2022
  delta_f_delta_t<- total_up - total_out
  
  
  # UPDATE ENVIRONMENT BLOCK
  # instantaneous concentrations of PFAS in fish tissue (at times step i)
  c_f<- c_f + (delta_f_delta_t*t_step)

  c_w0<-c_w*sin(seq_step[i]/10)
  if(c_w0 < 0){
    c_w0<-0  
  }
    
  # print(i)
  # print(c("total_up",total_up))
  # print(c("total_out",total_out))
  # print(c("delta_f_delta_t",delta_f_delta_t))
  # print(c("fish i ", results_c_f[i]))
  # print(c("kd", kd))
  # print(c("c_d", c_d))

  # append results of the time step to the vector 
  results_c_f<-append(results_c_f, c_f)
  results_c_w<-append(results_c_w, c_w0)

  results_c_s<-append(results_c_s, c_s)
  results_gill_up<-append(results_gill_up, gill_up)
  results_diet_up<-append(results_diet_up, diet_up)
  results_total_out<-append(results_total_out, total_out)

  if(is.na(c_f)){
    break
  }
  
}

length(results_c_w) == length(seq_step)
# results_c_f
data<-as.data.frame(results_c_f)
data$t<-c(0, seq_step)
data<-cbind(data, results_gill_up)
data<-cbind(data, results_diet_up)
data<-cbind(data, results_total_out)
data<-cbind(data, results_c_w)
data<-data[-1,]
  
data %>% 
  pivot_longer(cols = c(results_c_f, results_c_w)) %>% 
  ggplot(aes(t, value, color = name)) +
  geom_line()+
  facet_wrap(.~name, scales = "free_y")
  

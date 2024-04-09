# order is for c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")


inputFiles_list<-create_data_tables(
  species = c("Phy", "Pry", "Bgl", "Bas"), 
  group_species = c("plant", "fish", "fish", "fish"), 
  WB_kg = c(NA, 0.005, 0.145655, 0.495475), 
  m_O = c(1, 1, 1, 1),
  # diet = list(
  #     Pry = c(8.470000, 1115.000, 0.09800000, 0.15900000, 0,   0),
  #     Bgl = c(4.316667, 1844.167, 0.14233333, 0.14127500, 0.0, 0.0),
  #     Bas = c(2.705000, 2924.917, 0.07467473, 0.08748328, 0.0, 0.0)
  #     ), # order for PFAA
  diet = list(
      Pry = c(10.6, 1153,  0.194, 0.229,  0.029,   0),
      Bgl = c(35.9, 1408,  0.933, 0.704,  0.0754, 0.0),
      Bas = c(40.8, 3901,  0.125, 1.99,   0.0754, 0.0)
      ), # order for PFAA
  foodWeb = list(
      Phy = c(1.0, 0.0, 0.0, 0.0, 0.0), 
      Pry = c(0.05, 0.7, 0.0, 0.0, 0.0),
      Bgl = c(0.05, 0.3, 0.4, 0.0, 0.0),
      Bas = c(0, 0.4, 0.6, 0, 0.0)
      ),
  # foodWeb = list(
  #     Phy = c(1.0, 0.0, 0.0, 0.0, 0.0), 
  #     Pry = c(0.3, 0.7, 0.0, 0.0, 0.0),
  #     Bgl = c(0.3, 0.3, 0.4, 0.0, 0.0),
  #     Bas = c(0.2, 0.4, 0.6, 0, 0.0)
  #     ),
  C_WTO_ng_mL = c(0.1219425, 0.3862425, 0.03695603, 0.00372315, 0.0, 0.00), 
  # order is for c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
  C_s_ng_g = c(0.9623161, 0.2731140, 0.14765057, 11.04686791, 0.0, 0.00), # order of PFAA
  C_WTO_max_ng_mL = c(1.4930000, 0.8241000, 0.06314000, 0.00608600, 0.0, 0.00), # order of PFAA
  C_WTO_min_ng_mL = c(0.0166400, 0.1031000, 0.00971100, 0.00266000, 0.0, 0.00), # order of PFAA
  C_s_max_ng_g = c(3.5425000, 0.8300847, 0.45000000, 50.92500000, 0.0, 0.00), # order of PFAA
  C_s_min_ng_g = c(0.1534500, 0.1250000, 0.12500000,1.07676969, 0.0, 0.00),# order of PFAA
  C_OX = 8.0, 
  T = 14.0625,
  RMR = c(0, 0, 0, 0), # in mgO2/d whole fish
  min_diet = list(
      Pry = c(5.270, 705, 0.08060, 0.04030, 0,   0),
      Bgl = c(2.080, 380.5, 0.0396,  0.0396, 0.0, 0.0),
      Bas = c(0.03402, 2215, 0.03365, 0.03845, 0.0, 0.0)
      ), # order for PFAA
  max_diet = list(
      Pry = c(25.100, 1590, 0.91, 0.04030, 0,   0),
      Bgl = c(114.180, 4747, 3.833, 4.3830, 0.0, 0.0),
      Bas = c(63.0350, 4470, 0.16807, 1.6416, 0.0, 0.0)
)) # order for PFAA)


kRTable = read.csv('../Data/export_krTable.csv', row.names = "X", header = T)
colnames(kRTable)<-c("kr/kb", "PFAA", "chemID")

from run_polco import *
import time

TimoTimo("toy_model.1_split")
TimoTimo("toy_model.2_split")
TimoTimo("toy_model.3_split")
#TimoTimo("Ecoli_iAF1260_ox.glc.4_split")
#TimoTimo("iTM686.light.CLS.1_split")    


## this is optional, used for benchmarking
#~ f_out = open("running_times.txt","w")

#~ a = time.time()    
#~ for n in range(1,13):
    #~ TimoTimo("iTM686.glycogen.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("iTM686.glycogen",b-a))

#~ a = time.time()    
#~ for n in range(1,13):
    #~ TimoTimo("iTM686.light.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("iTM686.light",b-a))

#~ a = time.time()    
#~ for n in range(1,13):
    #~ TimoTimo("Ecoli_iAF1260_noox.glc.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("Ecoli_iAF1260_noox.glc",b-a))


#~ a = time.time()    
#~ for n in range(1,14):
    #~ TimoTimo("Ecoli_iAF1260_ox.glc.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("Ecoli_iAF1260_ox.glc",b-a))

#~ a = time.time()    
#~ for n in range(1,9):
    #~ TimoTimo("Lactis.glc.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("Lactis.glc",b-a))

#~ a = time.time()    
#~ for n in range(1,8):
    #~ TimoTimo("STherm.lcts.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("STherm.lcts",b-a))

#~ a = time.time()    
#~ for n in range(1,17):
    #~ TimoTimo("iNJ661.glyc.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("iNJ661.glyc",b-a))


#~ a = time.time()    
#~ for n in range(1,11):
    #~ TimoTimo("Ecoli_iJR904.mal_L.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("Ecoli_iJR904.mal_L",b-a))

#~ a = time.time()    
#~ for n in range(1,9):
    #~ TimoTimo("Ecoli_iJR904.glc.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("Ecoli_iJR904.glc",b-a))

#~ a = time.time()    
#~ for n in range(1,12):
    #~ TimoTimo("Ecoli_iJR904.fum.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("Ecoli_iJR904.fum",b-a))


#~ a = time.time()    
#~ for n in range(1,8):
    #~ TimoTimo("Ecoli_iJR904.trp_L.{0:d}_split".format(n))
#~ b = time.time()
#~ f_out.write("{0:s}\t{1}".format("Ecoli_iJR904.trp_L",b-a))

#~ f_out.close()





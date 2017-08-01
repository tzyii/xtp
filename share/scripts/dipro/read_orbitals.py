#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#



def readOrb(path,molIdx,nsaos_comp):

  molcoef=[]
  molcoef_add=[]
  coefs_full=[]
  LastLine=[]


  if molIdx!=1 and molIdx!=2:
     print "molIdx can only be 1 or 2!"
     sys.exit()
  

  molFile=open(path,'r')
  readMo='false'
  for line in molFile:
	  #skip comment lines
          if line.find('#') == 0:
	    continue

	  #stop reading when $end statement is reached
   	  if '$end' in line:
     	    break

	  #analyse mo-file header
	  if 'scfconv' in line:
	    cols=line.split('=')
	    str0=cols[1].split(' ')
	    scfconv=int(str0[0])
	  line = line.strip()

	  #read eigenvalue and calculate size of mo-block
	  if 'nsaos' in line:
	      lineCounter=1
	      readMo='true'
      	      cols=line.split('=')
	      nsaos=int(cols[2])
	      str1=cols[1].split(' ')
	      eigenvalue=[float(str1[0].replace('D','e'))]
	      ElementsInLastLine=nsaos%4
	      if ElementsInLastLine != 0:
 		  NumberOfLines=(nsaos/4+1)
	      else:
		  NumberOfLines=(nsaos/4)
	      continue

  	  #read mo-coefficients    
	  if readMo == 'true':     
  	    CoeffString=line
  	    #put the mo-coefficients into molcoef1 
  	    if lineCounter < NumberOfLines:
  	      for j in range(4):
  		molcoef.append(  CoeffString[0+j*20:20+j*20] )
  	      lineCounter+=1
  	    elif lineCounter == NumberOfLines and not 'nsaos' in line:
  	     #take care for non-complete lines
	     if ElementsInLastLine != 0:
  	       for k in range(ElementsInLastLine):
  	          molcoef.append(  CoeffString[0+k*20:20+k*20] )
	     else:	
  	       for k in range(4):
  	          molcoef.append(  CoeffString[0+k*20:20+k*20] )

	     for j in range(nsaos_comp):
  	     #generate field with zeros for the other molecule
  		molcoef_add.append( '0.00000000000000D-00' )
  	     #now glue eigenvalue, coefficients and zeros together
             if molIdx == 1:
  	        eigenvalue.extend(molcoef)	      
  	        eigenvalue.extend(molcoef_add)
	     else:
  	        eigenvalue.extend(molcoef_add)
  	        eigenvalue.extend(molcoef)	      
  	     #store complete mo into the mo-vector list	      
  	     coefs_full.append( eigenvalue )
  	     #re-initialize for next pass
  	     molcoef=[]
  	     molcoef_add=[]
  	     eigenvalue=[]
	     readMo='false'
  molFile.close()

  return coefs_full



def getOrbDim(path):
  """ extract "scfconv" and "nsaos" from the orbital read """
     


  molFile=open(path,'r')
  for line in molFile:
	  #skip comment lines
          if line.find('#') == 0:
	    continue

	  #stop reading when $end statement is reached
   	  if '$end' in line:
	    print "we should never come here!"
     	    break

	  #analyse mo-file header
	  if 'scfconv' in line:
	    cols=line.split('=')
	    str0=cols[1].split(' ')
	    scfconv=int(str0[0])

	  #read size of mo-block
	  if 'nsaos' in line:
	      lineCounter=1
	      readMo='true'
      	      cols=line.split('=')
	      nsaos=int(cols[2])
	      break

  molFile.close()

  return scfconv,nsaos	      
	      
	      
def writeMo(scfconv,nsaos,coefs_full,name):
  """coefs_full is the mo vector field obtained by readOrb, name should be alpha,beta or mos"""
 
  import sys
   
  outFile=open(name,'w')

  

  if name == "alpha":
    outFile.write("$uhfmo_alpha    scfconv=%d   format(4d20.14)\n#generated by merge_mos.py\n#\n" % (scfconv))
  if name == "beta":
    outFile.write("$uhfmo_beta    scfconv=%d   format(4d20.14)\n#generated by merge_mos.py\n#\n" % (scfconv))
  if name == "mos":
    outFile.write("$scfmo    scfconv=%d   format(4d20.14)\n#generated by merge_mos.py\n#\n" % (scfconv))

  ElementsInLastLineNew=nsaos % 4
  for i in range(nsaos):
  #loop over mos
     outFile.write("%6d  a      eigenvalue=%19.14lfD+00   nsaos=%d\n" % (i+1,coefs_full[i][0],nsaos))
     for j in range(nsaos/4):
     #loop over lines
  	outFile.write("%s%s%s%s\n" % (coefs_full[i][1+j*4],coefs_full[i][2+j*4],coefs_full[i][3+j*4],coefs_full[i][4+j*4]))
     if ElementsInLastLineNew > 0:
  	LastLine=[]
  	for k in range(ElementsInLastLineNew):
  	   #loop for elements in last line
  	   LastLine.append(coefs_full[i][k+1+(j+1)*4])
  	str3=''.join(LastLine)
  	outFile.write("%s\n" % (str3))
  outFile.write('$end\n')
  outFile.close()   
  
  
  
  
  
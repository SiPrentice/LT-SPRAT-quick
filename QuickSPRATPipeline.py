### Extracts spectrum from objects with reduced spectra from SPRATpipeline
### v.28-01-20 Includes the newest SPRATCorrections and acconts for object names that begin with letters
### v.18/3/20 now compares spectra scaled to the flux at a specific wavelength

'''
Extracts spectra from SPRAT pipeline and applies an airmass dependent correction to the flux
Required folders listed for input, output, and plots.
'''




import os
import glob
from pyraf import iraf
import pylab as plt
import numpy as np
from astropy.io import fits

###################

z=0.0
plotgal='y'
plotvel='n'
velocity=15500
E = 1e-5
Rv=3.1

input_location = '/InputSpectra/'
output_location='/OutputSpectra/'
plot_save_location = '/plots/'

corloc='./SPRATFluxCorrection.txt'


#####################

def dop(lamr,vel):
    vel=vel*-1
    a = (1.+vel/3.e5)**0.5
    b=(1.-vel/3.e5)**0.5
    return lamr*a/b

def dered(wav,flux,Rv,Ebv):
   lam=wav*0.0001
   Av=Ebv*Rv
   x=1/lam
   y=x-1.82
   a=1+(0.17699*y)-(0.50477*y**2)-(0.02427*y**3)+(0.72085*y**4)+(0.01979*y**5)-(0.77530*y**6)+(0.32999*y**7)
   b=(1.41338*y)+(2.28305*y**2)+(1.07233*y**3)-(5.38434*y**4)-(0.62251*y**5)+(5.30260*y**6)-(2.09002*y**7)
   AlAv=a+b/Rv
   Al=AlAv*Av
   #print Al
   F=10**(Al/2.5)
   delF= flux*F
   return delF

def name(fits_file):
    hdul = fits.open(fits_file)
    hdr = hdul[0].header
    if 'OBJECT' in hdr:

        return hdr['OBJECT'], hdr['DATE-OBS'], hdr['MJD'],hdr['AIRMASS']

        return hdr['OBJECT']

            
    else:
        print ('No defined object name in fits header -- please provide an object name\n') 
        return raw_input('Object name: ')

    
def header(fits_file):
    hdul = fits.open(fits_file)
    hdr = hdul[0].header
    return hdr

##################

dates=[]

files = glob.glob('.'+input_location+'*.fits')
for spectrum in files:
 justplot= False

 SN,dateobs,mjd,airmass = name(spectrum)[0],name(spectrum)[1],name(spectrum)[2],name(spectrum)[3]
 
 date=os.path.basename(spectrum)[4:12]
 if '_1_' in spectrum:
    	date = date+'_1'
 if '_2_' in spectrum:
    	date = date+'_2'
 if '_3_' in spectrum:
    	date = date+'_3' 
 print '\ndate = ', date
 
 

 file_list = glob.glob('.'+output_location+'*.txt')
 file_list.sort()
 if len(file_list)>0:
    for f in file_list:
        if date in f:
            justplot=True   

 if justplot == False:
    
    print 'Extracting spectrum'
    
    np.savetxt('.'+input_location+'allspec',[spectrum+'[4]'],fmt="%s")
        
    ## begin IRAFing

    ### Extract the 1D spectrum
    iraf.scopy(spectrum+'[4]', '.'+input_location+'allspeccal2')

    iraf.wspectext('.'+input_location+'allspeccal2', '.'+input_location+''+SN+'_'+date+'.txt',header='no')
    print 'Removing allspeccal2.fits'
    os.remove('.'+input_location+'allspeccal2.fits')
 else:
    print 'Skipping spectrum extraction'

 ### Correct the spectrum
 
 corrections= glob.glob('./CorrectionFiles/*.txt')
 clist = [ float(os.path.basename(c)[:-4])  for c in corrections]
 correction_index = np.argmin(abs(np.array(clist)-airmass))
 corloc = './CorrectionFiles/'+str(clist[correction_index])+'.txt'
 print 'Correcting flux for airmass', airmass, 'using', str(clist[correction_index])+'.txt'
    
 s = np.loadtxt(corloc,unpack=True, usecols=(0,1))
 txtlist = glob.glob('.'+input_location+'*.txt')
    
    
 for f in txtlist:
    if '.w.' not in f:
        if 'SPRAT' not in f:
            if date in f:
                file_to_use = f
    
 print 'Using', file_to_use
 a=np.loadtxt(file_to_use,unpack=True)
    
 for q in range(len(a[0])):


##################
justplot= False

file_list = glob.glob('./*.fits')
SN = name(file_list[0])

file_list = glob.glob('./*.txt')
file_list.sort()
if len(file_list)>0:
    for f in file_list:
        if SN in f:
            justplot=True

if justplot == False:
    
    print 'Extracting spectrum'
    
    ## create the fits list
    file_list = glob.glob('./*.fits');file_list.sort()
    final=[]
    for j in range(len(file_list)):
            final.append(file_list[j]+'[4]')
            final.append('')
    np.savetxt('./allspec',final,fmt="%s")

    ## extract the date
    date=file_list[0][6:14]
    
    print 'date = ', date
    
    ## begin IRAFing

    ### Extract the 1D spectrum
    iraf.scopy('@allspec', 'allspeccal2')

    iraf.wspectext('allspeccal2', SN+'_'+date+'.txt',header='no')
    print 'Removing allspeccal2.fits'
    os.remove('./allspeccal2.fits')
else:
    print 'Skipping spectrum extraction'

### Correct the spectrum    
s = np.loadtxt(corloc,unpack=True, usecols=(0,1))
txtlist = glob.glob('./*.txt')
    
    
for f in txtlist:
    if '.w.' not in f:
        if 'SPRAT' not in f:
            file_to_use = f
    
print 'Using', file_to_use
a=np.loadtxt(file_to_use,unpack=True)
    
for q in range(len(a[0])):

    find = np.argmin(abs(s[0]-a[0][q]))
    a[1][q] = a[1][q]/s[1][find]

         

 
 master = zip(a[0],a[1])
    
 ### Save the corrected spectrum
 np.savetxt('.'+output_location+os.path.basename(file_to_use)[:-4]+'.w.txt',master,fmt="%s", header=str(dateobs)+' '+ str(mjd))
 dates.append([os.path.basename(file_to_use)[:-4]+'.w.txt', dateobs,mjd])   
  
 ### Load the corrected spectrum
 n = np.loadtxt('.'+output_location+''+os.path.basename(file_to_use)[:-4]+'.w.txt',dtype='float',unpack=True)
    
 ### Plot the corrected spectrum dereddened for E
 m = np.max(a[1])
 a[1] = dered(a[0],a[1],Rv,E)
 plt.plot(a[0]/(1.+z),a[1]/max(a[1]),color='k',linewidth=1,label=SN+'\n$z=$'+str(z)+'\n$E_\mathrm{MW} = $'+str(E)[:4]+' mag')
    
 ### if necessary plot the galaxy lines
 if plotgal != 'n':

print SN, 'Corrected to SPRATFluxCorrections.txt'
master = zip(a[0],a[1])
    
### Save the corrected spectrum
np.savetxt(file_to_use[:-4]+'.w.txt',master,fmt="%s")
    
### Load the corrected spectrum
n = np.loadtxt(file_to_use[:-4]+'.w.txt',dtype='float',unpack=True)
    
### Plot the corrected spectrum dereddened for E
m = np.max(a[1])
a[1] = dered(a[0],a[1],Rv,E)
plt.plot(a[0]/(1.+z),a[1]/max(a[1]),color='k',linewidth=1,label=SN+'\n$z=$'+str(z)+'\n$E_\mathrm{MW} = $'+str(E)[:4]+' mag')
    
### if necessary plot the galaxy lines
if plotgal != 'n':

    l = [6548,6583, 4959,5007,5890,6717,6731]
    for line in l:
        ys = np.linspace(0,m*1.10)
        xs = [line for op in ys]
        plt.plot(xs,ys,color='grey',linestyle='dotted',zorder=0,linewidth=0.7)
    l = [6563,4861,4341,4100]
    for line in l:
        ys = np.linspace(0,m*1.10)
        xs = [line for op in ys]
        plt.plot(xs,ys,color='red',linestyle='dashed',zorder=0,linewidth=0.7)
        

 if plotvel !='n':

if plotvel !='n':

    l = [6355]
    for line in l:
        ys = np.linspace(0,m*1.10)
        xs = [dop(line,velocity) for op in ys]
        plt.plot(xs,ys,color='green',linestyle='dashed',zorder=0,linewidth=0.7)
                    
    

 plt.legend(loc='upper right')
 plt.ylim([0,1.1])
 plt.xlabel('Rest-frame wavelength [$\AA$]')
 plt.ylabel('Scaled flux')
 plt.savefig('.'+plot_save_location+SN+date+'.pdf',bbox_inches='tight')
 plt.close()
    
    
 ### Compare the two spectra
 a = np.loadtxt(file_to_use,unpack=True)
 scale = a[1][np.argmin(abs(a[0]-5000))]
 plt.plot(a[0]/(1.+z),a[1]/scale,color='k',alpha=1,linewidth=1,label='Original')

 a = np.loadtxt('.'+output_location+os.path.basename(file_to_use)[:-4]+'.w.txt',unpack=True)
 scale = a[1][np.argmin(abs(a[0]-5000))]
 plt.plot(a[0]/(1.+z),a[1]/scale,color='red',alpha=0.7,linewidth=1,label='Corrected')
 plt.legend()


 plt.xlabel('Rest-frame wavelength [$\AA$]')
 plt.ylabel('Scaled flux')
 plt.savefig('.'+plot_save_location+SN+date+'_compare.pdf',bbox_inches='tight')
 plt.close()
 
np.savetxt('./SpectraDates.txt',dates,fmt="%s")

plt.legend(loc='upper right')
plt.ylim([0,1.1])
plt.xlabel('Rest-frame wavelength [$\AA$]')
plt.ylabel('Scaled flux')
plt.savefig('./'+SN+'.pdf',bbox_inches='tight')
plt.close()
    
    
### Compare the two spectra
a = np.loadtxt(file_to_use,unpack=True)
plt.plot(a[0]/(1.+z),a[1]/max(a[1]),color='k',alpha=1,linewidth=1,label='Original')

a = np.loadtxt(file_to_use[:-4]+'.w.txt',unpack=True)
plt.plot(a[0]/(1.+z),a[1]/max(a[1]),color='red',alpha=0.7,linewidth=1,label='Corrected')
plt.legend()


plt.xlabel('Rest-frame wavelength [$\AA$]')
plt.ylabel('Scaled flux')
plt.savefig('./'+SN+'_compare.pdf',bbox_inches='tight')
plt.close()


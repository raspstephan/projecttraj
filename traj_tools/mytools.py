import sys, struct
from numpy import *
from subprocess import call, Popen, PIPE
from datetime import datetime, timedelta
import warnings
#import Nio
from scipy import weave
import os

def seq(start, end, step=1.):
  return linspace(start, end, round((end-start)/step)+1)
  
def ran(start, end, step):
  return range(start,end+1,step)
  
def crosscorr(ZR1,ZR2,maxlag):
  #Berechnet die Crosscorrel zw. den gleich langen Zeitreihen 1 und 2
  #Es werden [-maxlag,...,maxlag] Wert berechnet
  #Posit. lag: ZR1 fuehrt
  n=len(ZR1)
  if len(ZR2)!=n:
    print('error')
    
  #bei neg. dt: ZR1 vorne kuerzen, ZR2 hinten
  dtint=[]; erg=[]
  for it in xrange(2*maxlag+1):
    dt=it-maxlag
    dtint.append(dt)
    if dt>=0:
      erg.append(corrcoef(ZR1[:n-dt],ZR2[dt:])[0,1])
    else:
      erg.append(corrcoef(ZR1[-dt:],ZR2[:n+dt])[0,1])
  return array(erg),array(dtint)

def ZeitEinsetzen(Str,t):
  if isinstance(t,datetime):
    Str=Str.replace('%yyyy',str(t.date().year))
    Str=Str.replace('%mm','%0.2d' % (t.date().month))
    Str=Str.replace('%dd','%0.2d' % (t.date().day))
    Str=Str.replace('%hh','%0.2d' % (t.time().hour))
    Str=Str.replace('%ii','%0.2d' % (t.time().minute))
    Str=Str.replace('%ss','%0.2d' % (t.time().second))
  elif isinstance(t,timedelta):
    stot=int(t.total_seconds())
    d=stot/86400
    stot=stot-d*86400
    h=stot/3600
    stot=stot-h*3600
    i=stot/60
    stot=stot-i*60
    s=stot
    Str=Str.replace('%dd','%0.2d' % (d))
    Str=Str.replace('%hh','%0.2d' % (h))
    Str=Str.replace('%ii','%0.2d' % (i))
    Str=Str.replace('%ss','%0.2d' % (s))
  else:
    print('Error!')
  return Str
  
def SucheFileName(Maske,Stellen=4,dummyfile=True):
  #Sucht mit Hilfe der Maske den Dateinamen mit der kleinsten Nummer, der noch nicht existiert. Die Nummer muss in der Maske mit %num angegeben sein.
  #Optional kann die Anzahl der Stellen angegeben werden. Bei dummyfile=True wird sofort ein dummyfile mit dem Namen angelegt, um zu verhindern,
  #dass ein anderer Prozess den gleichen Namen findet.
  print('Warnung: Methode unzuverlaessig! Doppeltes Fileanlegen moeglich.')
  n=0
  fmtstr='%0xi'
  fmtstr=fmtstr.replace('x',str(Stellen))
  while True:
    try:
      numstr=fmtstr%(n)
      fid = open(Maske.replace('%num',numstr), 'r')
      fid.close()
    except:
      if dummyfile:
	fid = open(Maske.replace('%num',numstr), 'w')
	fid.close()
      break
    n+=1
  return Maske.replace('%num',numstr)
  
def Zeitschritte(Start,Ende,Step):
  #Berechnet die Zeitschritte bei Schrittweite Step zw Start und Ende einschliesslich
  return int((Ende-Start).total_seconds()/Step.total_seconds())+1

#def LeseNioGrib(Filename,Varname,ilev=-1):
  ##Liest Gribfile mit Hilfe der PyNio-Bibliothek ein
  #grbfile=Nio.open_file(Filename+".grb1","r")
  #var=grbfile.variables[Varname]
  #if var.rank==3:
    #nlev=var.shape[0]
    #if nlev>1 and ilev==-1:
      #print("Error: Multilevel-Field. You have to specify the level")
  #else:
    #nlev=1
  #if ilev==-1:
    #ilev=0
  #if ilev>=nlev:
    #print("Error: Requested level not present in Gribfile.")
  
  #if var.rank==3:
    #Feld=array(var[ilev].transpose(),float)
  #else:
    #Feld=array(var[:].transpose(),float)
  #nlon=Feld.shape[0]; nlat=Feld.shape[1]
  
  #grbfile.close()
  
  #return Feld

def FindeRecNr(GribFile,DefStr):
  #Behandelt unsortierte Gribfiles ab COSMO Version 4.24. Es sucht im wgrib -v output nach
  #DefStr und gibt die entsprechende RecNr zurueck
  
  wgribcommand='wgrib -v '+GribFile+' | grep "%s"'%(DefStr)
  
  p=Popen(wgribcommand,shell=True,stdout=PIPE)
  output, errors = p.communicate()
  
  #Mehrere Recs gefunden
  if(output.count('\n')>1):
    raise Exception('Mehrere passende Records gefunden')
  if(output==''):
    raise Exception('Kein passender Record gefunden')
  
  RecNr=int(output[0:output.find(':')])
  
  return RecNr

def LeseGrib(Filename,RecNr,nlon,nlat):
  
  #RecNr kann auch eine String enthalten der die Variable mit obiger Funktion identifiziert.
  if (type(RecNr) is str):
    RecNr=FindeRecNr(Filename,RecNr)
  
  #tmpfile generieren
  compname=os.uname()[1]
  pID=str(os.getpid())
  tmpfile='/local/tobias.selz/tmpfile_'+compname+'_'+pID
  
  try:
    nlev=len(RecNr)
  except:
    nlev=1
    RecNr=[RecNr]
    
  if nlev==1:
    Feld=zeros((nlon,nlat),'float')
  else:
    Feld=zeros((nlon,nlat,nlev),'float')
  
  for ilev in xrange(nlev):
    call('rm -f '+tmpfile,shell=True)
    wgribcomand='wgrib -d %i %s -o %s' % (RecNr[ilev],Filename,tmpfile)
    p=Popen(wgribcomand,shell=True,stdout=PIPE)
    output, errors = p.communicate()
    #File laden
    f=open(tmpfile,'r+b')
    dummy=f.read(4)
    for ilat in xrange(nlat):
      dummy=f.read(4*nlon)
      if nlev==1:
        Feld[:,ilat]=struct.unpack(str(nlon)+'f',dummy)
      else:
        Feld[:,ilat,ilev]=struct.unpack(str(nlon)+'f',dummy)
    dummy=f.read(4)
    f.close()
    call('rm -f '+tmpfile,shell=True)
    
  return Feld
  
def ErsetzeGribJahr(Infile,Outfile,neuesJahr):
  #Ersetzt das Jahr im GribRec durch neuesJahr. Beachte: Jahr=Jahr des Jahrhunderts
  #Anzahl der prozessierten Recs wird zurueckgegeben
  
  #File laden
  f=open(Infile,'r+b')
  g=open(Outfile,'w+b')
  RecAnz=0
  while(True):
    #Grib-Sektion 0 lesen
    Str=f.read(8)
    if len(Str)<4:
      break #EOF
    if Str[0:4]!='GRIB':
      print('Error: kein GRIB!')
    RecAnz+=1
    #Recordlaenge bestimmen
    RecLen=struct.unpack('>I', '\x00'+Str[4:7])
    #Ganzen Record lesen
    f.seek(-8,1)
    clist=list(f.read(RecLen[0]))
    #Jahr ersetzen
    clist[20]=struct.pack('B',neuesJahr)
    #Record schreiben
    g.write(''.join(clist))
  g.close()
  f.close()
  return RecAnz
    
def SchreibeDat(Filename,Feld,res='f',mode='w'):
  #Schreibt Feld und Filename auf die Platte. Res gibt die Aufloesung an: f=4Byte, d=8Byte.
  #Mode gibt an, ob bestehende Files ueberschrieben werden sollen (w) oder angehaengt (a)
  dim=len(Feld.shape)
  if dim>3:
    raise Exception('dim muss <=3 sein')
  
  #Resolution pruefen, BytesPerValue setzen
  if res!='f' and res!='d':
    raise Exception('res muss f oder d sein')
  
  nlon=Feld.shape[0]
  
  if dim==1:
    nlat=1; nlev=1
  elif dim==2:
    nlat=Feld.shape[1]; nlev=1
  else:
    nlat=Feld.shape[1]; nlev=Feld.shape[2]
  
  f=open(Filename,mode);
  for ilev in xrange(nlev):
    for ilat in xrange(nlat):
      for ilon in xrange(nlon):
        if dim==1:
          f.write(struct.pack(res,Feld[ilon]))
        elif dim==2:
          f.write(struct.pack(res,Feld[ilon,ilat]))
        else:
          f.write(struct.pack(res,Feld[ilon,ilat,ilev]))
  f.close()
  return 0

def LeseDat(Filename,Feld,skip=0,res='f'):
  #Fuellt das vorgegebene Feld linear mit Daten aus dem File auf
  #Reihenfolge:lons,lats,levs
  #Feld kann 1-3 dimensional sein
  #Bei Angabe von skip wird skip mal die Feldgroesse uebersprungen
  #Ist dim(Feld)=3 darf skip nicht angegeben werden
  #res: Resolution auf der Platte: f=float=4Byte, d=double=8Byte
  
  #Resolution pruefen, BytesPerValue setzen
  if res!='f' and res!='d':
    raise Exception('res muss f oder d sein')
  if res=='f':
    BpV=4
  else:
    BpV=8
    
  dim=len(Feld.shape)
  if dim>3:
    raise Exception('Feld darf hoechstens 3d sein')
    
  nlon=Feld.shape[0]
  if dim==1:
    nlat=1; nlev=1
  elif dim==2:
    nlat=Feld.shape[1]; nlev=1
  else:
    nlat=Feld.shape[1]; nlev=Feld.shape[2]
  #if skip!=0 and dim==3:
    #raise Exception('Skip wird bei 3d-Feldern nicht unterstuetzt')
    
  f=open(Filename,'r+b')
  f.seek(skip*nlon*nlat*nlev*BpV)
  for ilev in xrange(nlev):
    for ilat in xrange(nlat):
      dummy=f.read(BpV*nlon)
      if dim==1:
        Feld[:]=struct.unpack(str(nlon)+res,dummy)
      if dim==2:
        Feld[:,ilat]=struct.unpack(str(nlon)+res,dummy)
      if dim==3:
        Feld[:,ilat,ilev]=struct.unpack(str(nlon)+res,dummy)
  f.close()
  return 0
  
def intC2AGrid(Feld,var):
  #Interpoliert die Windvariablen u und v vom C auf ein A Gitter. Die Dimension des Felds bleibt unveraendert,
  #d.h. es werden die westlichsten und suedlichsten Werte kopiert.
  if var=='u':
    Feld[1:,:]=(Feld[:-1,:]+Feld[1:,:])/2.
  elif var=='v':
    Feld[:,1:]=(Feld[:,:-1]+Feld[:,1:])/2.
  else:
    raise Exception('var muss u oder v sein.')
  return 0
  
  
#def LeseDat(Filename,ilev,nlon,nlat):
  
  #if ilev!=0:
    #print('Error! Noch nicht implementiert!')
  #Feld=zeros((nlon,nlat),float)
  #f=open(Filename,'r+b')
  #for ilat in xrange(nlat):
    #dummy=f.read(4*nlon)
    #Feld[:,ilat]=struct.unpack(str(nlon)+'f',dummy)
  #f.close()
  #return Feld

def BerechneCosmoRefState(z,psl=100000.,Tsl=288.15,beta=42.):
  #Berechnet p0 , rho0 und T0 auf Voll-Modellflaechen der COSMO Base State Atmosphaere gemaess Gl 3.85 der Doku. Die Standardparameter sind bereits vorgegeben
  #Uebergeben werden muss die Hoehe der Halbflaechen z
  
  #z checken, Dimension ermitteln
  if len(shape(z))!=3:
    raise Exception('z muss 3d sein!')
  #Konstanten
  Rd=287.05; g=9.80665
  #p0 auf Halblevels berechnen (4.11)
  p0hl=psl*exp(-Tsl*(1.-sqrt(1.-2*beta*g*z/(Rd*Tsl**2)))/beta)
  #rho0 auf Vollflaechen berechnen (4.12)
  rho0=1./g*(p0hl[:,:,1:]-p0hl[:,:,:-1])/(z[:,:,:-1]-z[:,:,1:])
  #p0 auf Vollflaechen interpolieren (4.14)
  p0=0.5*(p0hl[:,:,1:]+p0hl[:,:,:-1])
  #T0 auf Vollflaechen berechnen (4.16)
  T0=p0/Rd/rho0
 
  return p0,rho0,T0

def BerechneCrosskorr_fft(f,g):
  #Berechnet die Crosscorrelationen der Felder f und g mittels FFT. f und g koennen 1- oder 2dim sein und muessen gleich gross sein.
  #Berechnet wird die inverse FFT des Produktes der Fouriertransformierten von f und g, letztere komplex konjugiert (Wiener-Chintschin-Theorem)
  #Es wird nicht normiert oder Mittelwert entfernt. Das muss man vorher selbst machen.
  sf=f.shape
  sg=g.shape
  dim=len(sf)
  if dim<1 or dim>2:
    raise Exception('Felder muessen 1- oder 2-dimensional sein.')
  #**********************************
  #1d-Version
  if dim==1:
    if sf[0]!=sg[0]:
      raise Exception('Beide Felder muessen die gleiche Dimension haben.')
    #Feldgroesse  
    Nx=sf[0]
    #CrossKorrelation
    M=real(fft.fftshift(fft.ifft(fft.fft(f)*conj(fft.fft(g)))))/Nx
  #***********************************  
  #2d-Version
  if dim==2:
    if sf[0]!=sg[0] or sf[1]!=sg[1]:
      raise Exception('Beide Felder muessen die gleiche Dimension haben.')
    #Feldgroesse  
    Nx=sf[0]; Ny=sf[1]
    #CrossKorrelation
    M=real(fft.fftshift(fft.ifft2(fft.fft2(f)*conj(fft.fft2(g)))))/Nx/Ny
    
  return M
  
def BerechneCrosskorr_shift(f,g,maxshift=0):
  #Berechnet alle moeglichen Crosscorrelationen der Felder f und g mittels Verschiebung. f und g koennen 1- oder 2dim sein und muessen gleich gross sein.
  #Berechnet wird M[ix,iy]=mean(f[x,y]*g[x+ix,y+iy]), positives ix,iy verschiebt g gegenueber f nach rechts und nach oben.
  #Es wird nicht normiert oder Mittelwert entfernt. Das muss man vorher selbst machen.
  #Maxshift gibt an um wie viele Gitterpunkte maximal verschoben werden soll. 0 steht fuer keine Grenze.
  sf=f.shape
  sg=g.shape
  dim=len(sf)
  if dim<1 or dim>2:
    raise Exception('Felder muessen 1- oder 2-dimensional sein.')
  #**********************************
  #1d-Version
  if dim==1:
    if sf[0]!=sg[0]:
      raise Exception('Beide Felder muessen die gleiche Dimension haben.')
    #Feldgroesse  
    Nx=sf[0];
    #Maxshift bestimmen
    if(maxshift==0):
      msx=Nx-1
    else:
      msx=maxshift
    #Autokorrmatrix anlegen
    M=zeros(2*msx+1,float)
    #Cpp code
    code = r"""
    for(int ix=0;ix<2*msx+1;ix++)
    {
      int dx=ix-msx;
      //Schleifenlaenge
      int nlx=Nx-abs(dx);
      //offsets berechnen
      int ofx=0, ogx=0;
      if(dx>=0) ofx=dx;
      else ogx=-dx;
      //Mittelungsschleife
      for(int i=0;i<nlx;i++) M(ix)+=f(i+ofx)*g(i+ogx);
      M(ix)/=nlx;
    }
    """
    err=weave.inline(code,['Nx','msx','f','g','M'],type_converters=weave.converters.blitz,compiler='gcc',headers=["<math.h>"])
  #***********************************  
  #2d-Version
  if dim==2:
    if sf[0]!=sg[0] or sf[1]!=sg[1]:
      raise Exception('Beide Felder muessen die gleiche Dimension haben.')
    #Feldgroesse  
    Nx=sf[0]; Ny=sf[1]
    #Maxshift bestimmen
    if(maxshift==0):
      msx=Nx-1; msy=Ny-1
    else:
      msx=int(maxshift[0]); msy=int(maxshift[1])
    #Autokorrmatrix anlegen
    M=zeros((2*msx+1,2*msy+1),float)
    #Cpp code
    code = r"""
    for(int ix=0;ix<2*msx+1;ix++)
    {
      int dx=ix-msx;
      for(int iy=0;iy<2*msy+1;iy++)
      {
	int dy=iy-msy;
	//Schleifenlaengen
	int nlx=Nx-abs(dx);
	int nly=Ny-abs(dy);
	//offsets berechnen
	int ofx=0, ofy=0;
	int ogx=0, ogy=0;
	if(dx>=0) ofx=dx;
	else ogx=-dx;
	if(dy>=0) ofy=dy;
	else ogy=-dy;
	//Mittelungsschleife
	for(int i=0;i<nlx;i++)
	{
	  for(int j=0;j<nly;j++)
	  {
	    M(ix,iy)+=f(i+ofx,j+ofy)*g(i+ogx,j+ogy);
	  }
	}
	M(ix,iy)/=nlx*nly;
      }
    }
    """
    err=weave.inline(code,['Nx','Ny','msx','msy','f','g','M'],type_converters=weave.converters.blitz,compiler='gcc',headers=["<math.h>"])
  
  return M
  
def Int_2DCorr(Corr2D,dx):
  #Berechnung die 1d-Korrfkt aus der zweidim. durch Ausintegrieren der Richtungsabhaengigkeit
  
  #dx und dy extrahieren
  if(len(dx)==1):
    dy=dx
  elif(len(dx)==2):
    dy=dx[1]; dx=dx[0]
  else:
    raise Exception('dx muss 1 oder 2dim sein.')
  
  #dr-Diskretisierungsintervall
  dr=max(dx,dy)
  Nx=Corr2D.shape[0]; Ny=Corr2D.shape[1]
  
  #Position des Ursprung
  i0=Nx/2; j0=Ny/2
  #Matrix mit Quadratabstaenden vom Ursprung erstellen
  ix=array(xrange(Nx))-i0; iy=array(xrange(Ny))-j0
  absMatrix=outer(ix*dx,ones(Ny,float))**2+outer(ones(Nx,float),iy*dy)**2

  #Anzahl der dr-Schritte bestimmen
  rmaxx=((Nx-1)/2+0.5)*dx; rmaxy=((Ny-1)/2+0.5)*dy
  rmax=min(rmaxx,rmaxy)
  N=int(rmax/dr-0.5)+1
  
  Corr1D=zeros(N,float)
  Corr1D[0]=Corr2D[i0,j0]
  #Schleife ueber alle k>0
  for ir in xrange(1,N):
    #Ringmaske fuer aktuelle ik berechnen
    Maske=(absMatrix>((ir-0.5)*dr)**2)&(absMatrix<=((ir+0.5)*dr)**2)
    #Eintraege in den Ringen mitteln
    Corr1D[ir]=(Corr2D*Maske).sum()/len(Maske.nonzero()[0])
  
  return [Corr1D,dr]
  
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#WICHTIG!!! Namen und Funktionsweise der Spektraltools muss ueberarbeitet werden!!!
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def DFT_1d(f,dx): #getestet! ok!
  #Berechnet f(k)=1/sqrt(L)*int_0^L dx f(x) exp(-ikx)
  #Zuvor wird ein Detrending durchgefuehrt. Der rechteste Punkt wird dadurch dem ersten gleichgesetzt.
  #Bei der Fouriertransformation bleibt er unberuecksichtigt. Diese wird also nur ueber N-1 Gitterpunkte 
  #ausgefuehrt.
  #Redundante Komponenten (Koeffizienten fuer grosse k bzw. neg. k) werden abgeschnitten.
  
  #Grundgroessen ermitteln
  N=len(f)
  L=float(dx*(N-1))
  dk=2*pi/L
  
  #k-Vektor erzeugen
  k=dk*array(xrange(N/2))
  
  #Trend entfernen, dazu f kopieren, damit das Original erhalten bleibt
  df=(f[-1]-f[0])/(N-1)
  fp=f.copy()
  fp-=df*xrange(N)
  
  #Transformieren, diskrete Form entsprechend renormieren.
  spf=sqrt(L)/float(N-1)*fft.fft(fp,N-1)
 
  return [spf[:len(k)],k]

  
def DFT_2d(f,dx):
  #Berechnet f(k)=1/(L_x*L_y)*int_0^L_x dx int_0^L_y dy f(x,y) exp(-ikx-ily)
  #Zuvor wird ein 2d-Detrending durchgefuehrt. Die rechtesten und obersten Punkt wird dadurch den linkesten und untersten gleichgesetzt.
  #Bei der Fouriertransformation bleiben die rechtesten und obersten unberuecksichtigt. Diese wird also nur ueber Nx-1,Ny-1 Gitterpunkte 
  #ausgefuehrt.
  #Die Koeffizientenmatrix wird verschoben, so dass kx,ky=0,0 in der Mitte liegt. Redundante Komponeneten werden beibehalten. Der Betrag
  #der Koeffizienten liegt punksymmetrisch um den Ursprung. 
  
  #Grundgroessen ermitteln
  s=f.shape
  if len(s)!=2:
    raise Exception('f muss zweidimansional sein')
  if len(dx)!=2:
    raise Exception('dx muss zweidimansional sein')
  Nx=s[0]; Ny=s[1]
  Lx=float(dx[0]*(Nx-1)); Ly=float(dx[1]*(Ny-1))
  dkx=2*pi/Lx; dky=2*pi/Ly
  
  #k-Vektoren erzeugen. Da das Spektrum spaeter verschoben wird muessen auch neg. k-Vektoren erzeugt werden.
  #Beachte: Nur Nx-1,Ny-1 Komponenten werden transformiert
  kx=dkx*array(xrange(-((Nx-1)/2),Nx/2)); ky=dky*array(xrange(-((Ny-1)/2),Ny/2))
  
  #Trend entfernen, f kopieren, damit Original erhalten bleibt
  fp=f.copy()
  #zunaechst x-Trends
  for j in xrange(Ny):
    df=(fp[-1,j]-fp[0,j])/float(Nx-1)
    fp[:,j]-=df*xrange(Nx)
  #y-Trends
  for i in xrange(Nx):
    df=(fp[i,-1]-fp[i,0])/float(Ny-1)
    fp[i,:]-=df*xrange(Ny)
    
  #Transformieren, diskrete Form entsprechend renormieren.
  spf=1.0/float((Nx-1)*(Ny-1))*fft.fft2(fp,[Nx-1,Ny-1])
  #spf=fft.fft2(fp,[Nx-1,Ny-1])
  spf=fft.fftshift(spf)
  
  return [spf,kx,ky]
  
def Int_2DSp(aspf,kx,ky):
  #Integriert ein 2d Spektrum und berechnet so ein eindimensionales. Berechnet wird:
  #E(k)=1/dk*sum aspf mit k im dk-Ring um k. Berechnung gemaess Whaite and Snyder.
  #Bei verschiedenen dkx und dky wird das Maximum als dk verwendet (wie bei Whaite, anders als bei Bierdel).
  #Ein groesseres dk fuehrt zu einer besseren Mittelung.
  #Es wird angenommen, dass aspf breits das Betragsquadrat der Spektralkoeffizienten enthaelt.
  #Im Fall der KE ist aspf=(abs(spu)**2+abs(spv)**2)/2
  #Der Ursprung muss in die Mitte der Domain verschoben sein.
  
  #dk und Nxy bestimmen
  dkx=kx[1]-kx[0]; dky=ky[1]-ky[0]
  dk=max(dkx,dky)
  Nx=len(kx); Ny=len(ky)
  
  #k-Vektor des Output-Spektrums bestimmen
  rmaxx=((Nx-1)/2+0.5)*dkx; rmaxy=((Ny-1)/2+0.5)*dky
  rmax=min(rmaxx,rmaxy)
  N=int(rmax/dk-0.5)+1
  k=dk*array(xrange(N))
 
  #Position des Ursprung
  ux=Nx/2; uy=Ny/2
  if(kx[ux]!=0 or ky[uy]!=0):
    raise Exception('Ursprung nicht zentriert.')
  #Matrix mit k-Quadratabstaenden vom Ursprung erstellen
  ix=array(xrange(Nx))-ux; iy=array(xrange(Ny))-uy
  absMatrix=outer(ix*dkx,ones(Ny,float))**2+outer(ones(Nx,float),iy*dky)**2

  #Output-Spektrum anlegen
  spf1d=zeros(N,float)
  spf1d[0]=1./dk*aspf[ux,uy]
  #Schleife ueber alle k>0
  for ik in xrange(1,N):
    #Ringmaske fuer aktuelle ik berechnen
    M=(absMatrix>(k[ik]-dk/2.)**2)&(absMatrix<=(k[ik]+dk/2.)**2)
    #Eintraege in den Ringen summieren
    spf1d[ik]=1./dk*(M*aspf).sum()
    #Schleife beenden, wenn der Ring den Rand erreicht hat
    #plt.contourf(M)
    #plt.show()
   
  return [spf1d,k]
   
def Spektrum1D(psi,dx):
  #Berechnet |psi(k)|^2 mit psi(k)=1/sqrt(L)*int_0^L dx psi(x) exp(ikx)
  #Zuvor wird ein Detrending durchgefuehrt
  
  #Grundgroessen ermitteln
  N=len(psi)
  L=float(dx*(N-1))
  dk=2*pi/L
  q=N/2+1 #Groesse des k-Vektors einschliesslich 0)
  kmax=(q-1)*dk
  
  #k- und Lambda-Vektor erzeugen
  k=dk*array(xrange(q+1))
  ld=2*pi/k
  
  #Trend entfernen
  dpsi=(psi[-1]-psi[0])/(N-1)
  psi-=dpsi*xrange(N)
  
  #Transformieren, diskrete Form entsprechend normieren.
  sp=float(N)/L*fft.fft(psi)
  
  #Betragsquadrat und redundanten Teil entfernen
  sp=L*abs(sp[:q+1])**2
  
  return [sp,k,ld]
  
  
  
def PowerDFT(Feld,dalpha):
  #Berechnet das Powerspektrum des 2D-Felds mit gleichfoermigem Gitterabstand dalpha in Grad
  #Dazu wird die diskrete Fouriertransformation verwendet
  
  #Dimension ermitteln
  s=Feld.shape
  if len(s)!=2:
    raise TypeError, 'Feld muss zweidimansional sein'
  nx=s[0]; ny=s[1]
  
  #dkxy ermitteln (fuer dkx ein mittleres dk)
  r_e=6371229.
  dkx=2*pi/(r_e*(nx-1)*deg2rad(dalpha)*cos((ny-1)*deg2rad(dalpha)/4.0))
  dky=2*pi/(r_e*(ny-1)*deg2rad(dalpha))
  #dk (k-Schritt des Spektrums) ist Minimum von beiden
  dk=min(dkx,dky)
  #nk ermitteln: Wie viele k-Werte bis zum Rand (DC-Komp nicht mitgezaehlt) koennen analysiert werden
  nkx=(nx-1)//2; nky=(ny-1)//2
  #max k-Radius ermitteln, bis zu dem integriert werden kann
  rk=max(nkx*dkx,nky*dky)
  #entsprechendes nk
  nk=int(rk/dk)
  
  #Fourier-Transformation
  f=abs(fft.fftshift(fft.fft2(Feld)))**2
  #Index der DC-Komponente
  nxdc=nx//2; nydc=ny//2
  
  #Ausgabefelder initialisieren
  power=zeros(nk+1,float) #+1 fuer DC-Komp
  npower=zeros(nk+1,int)
  
  #Punkte durchgehen (obere Haelfte)
  for iy in xrange(nydc,ny):
    for ix in xrange(nx):
      k=sqrt(((ix-nxdc)*dkx)**2+((iy-nydc)*dky)**2)
      idx=int(k/dk+0.5)
      if idx>nk:
	continue
      if iy==nydc:
	power[idx]+=f[ix,iy]
	npower[idx]+=1
      else:
	power[idx]+=2*f[ix,iy]
	npower[idx]+=2
  #Mittelwerte bilden
  power*=1.0/npower
  for ik in xrange(1,nk+1):
    power[ik]*=dk*ik*2*pi
    
  #k und lambda-Vektoren erzeugen
  k=arange(0.,nk+0.5,1)*dk
  lambd=2*pi/k
  
  return [power,k,lambd]
  
def PowerDCT(Feld,dalpha):
  #Berechnet das Powerspektrum des 2D-Felds mit gleichfoermigem Gitterabstand dalpha in Grad
  #Dazu wird die diskrete Cosinustransformation verwendet
  
  #Dimension ermitteln
  s=Feld.shape
  if len(s)!=2:
    raise TypeError, 'Feld muss zweidimansional sein'
  nx=s[0]; ny=s[1]
  
  #dkxy ermitteln (fuer dkx ein mittleres dk)
  r_e=6371229.
  dkx=pi/(r_e*(nx-1)*deg2rad(dalpha)*cos((ny-1)*deg2rad(dalpha)/4.0))
  dky=pi/(r_e*(ny-1)*deg2rad(dalpha))
  
  #dk (k-Schritt des Spektrums) ist Minimum von beiden
  dk=min(dkx,dky)
  #nk ermitteln: Wie viele k-Werte bis zum Rand (DC-Komp nicht mitgezaehlt) koennen analysiert werden
  nkx=nx-1; nky=ny-1
  #max k-Radius ermitteln, bis zu dem integriert werden kann
  rk=max(nkx*dkx,nky*dky)
  #entsprechendes nk
  nk=int(rk/dk)
  
  #Cosinus-Transformation
  f=dct2(Feld)**2
  
  #Ausgabefelder initialisieren
  power=zeros(nk+1,float) #+1 fuer DC-Komp
  npower=zeros(nk+1,int)
  
  #Punkte durchgehen
  #for iy in xrange(ny):
    #for ix in xrange(nx):
      #k=sqrt((ix*dkx)**2+(iy*dky)**2)
      #idx=int(k/dk+0.5)
      #if idx>nk:
	#continue
      #power[idx]+=f[ix,iy]
      #npower[idx]+=1
  ##Mittelwerte bilden
  #power*=1.0/npower
  #for ik in xrange(1,nk+1):
    #power[ik]*=dk*ik*2*pi
    
  #C-Routine
  code = r"""
  for(int iy=0; iy<ny; iy++)
  {
    for(int ix=0; ix<nx; ix++)
    {
      double k=sqrt(pow(double(ix)*double(dkx),2)+pow(double(iy)*double(dky),2));
      int idx=int(k/double(dk)+0.5);
      if(idx>nk) continue;
      power(idx)+=f(ix,iy);
      npower(idx)++;
    }
  }
  for(int ik=1;ik<=nk;ik++) power(ik)*=double(dk)*double(ik)*2.0*3.1415927/npower(ik);
  power(0)/=npower(0);
  """
  err=weave.inline(code,['nx','ny','dkx','dky','dk','nk','f','power','npower'],type_converters=weave.converters.blitz,compiler='gcc',headers=["<math.h>"])
  
  #k und lambda-Vektoren erzeugen
  k=arange(0.,nk+0.5,1)*dk
  lambd=2*pi/k
  
  return [power,k,lambd]
  
def PowerSp2D(Feld,dx):

  #Dimension ermitteln
  s=Feld.shape
  if len(s)!=2:
    raise TypeError, 'Feld muss zweidimansional sein'
  nx=s[0]; ny=s[1]

  #dk ermitteln
  dkx=2*pi/(nx*dx)
  dky=2*pi/(ny*dx)
  dk=min(dkx,dky)
  
  #idxmax ermitteln: kmax=idxmax*dk
  idxmax=(min(nx,ny)-1)//2
  
  #Nullen initialisieren
  power=zeros(idxmax+1,float)
  npower=zeros(idxmax+1)

  #Fourier-Transformation
  kFeld=fft.fft2(Feld)
  #Betragsquadrat bilden
  kFeld=abs(kFeld)**2
  
  #Loop in C++-Code
  code = r"""
  double k;
  for(int ix=0;ix<nx;ix++)
  {
    for(int iy=0;iy<ny;iy++)
    {
      k=sqrt(pow(ix*dkx,2)+pow(iy*dky,2));
      int idx=int(k/dk+0.5);
      if(idx>idxmax) continue;
      power(idx)+=kFeld(ix,iy);
      npower(idx)+=1;
    }
  }
  """
  err=weave.inline(code,['dkx','dky','dk','nx','ny','idxmax','kFeld','power','npower'],type_converters=weave.converters.blitz,compiler='gcc',headers=["<math.h>"])
  
  #Loop in Python-Code
  #for ix in xrange(nx):
    #for iy in xrange(ny):
      ##Skalares k berechnen
      #k=sqrt((ix*dkx)**2+(iy*dky)**2)
      ##Zu k gehoerenden Index berechnen: [0,0.5dk[->0, [0.5dk, 1.5dk[->1,...
      #idx=round(k/dk);
      #if idx>idxmax:
        #continue
      #power[idx]+=pFeld[ix,iy]
      #npower[idx]+=1

  #k und lambda-Vektoren erzeugen
  k=arange(0.,idxmax+0.5,1)*dk
  lambd=2*pi/k
  
  #Multiplizieren mit Integrationsflaeche, div durch Anzahl
  #k*fft(...) wegen Powerspektrum
  power*=2*pi*k/npower

  return [power, k, lambd]

def SpFilterDCT(Feld,dalpha,Art,Int):

  #Eingabe ueberpruefen
  if Art!='l' and Art!='k':
    raise ValueError, 'Art muss k oder l sein'
  if len(Int)!=2:
    raise TypeError, 'Int muss 2dim Vektor sein'
  if Int[0]>Int[1]:
    raise ValueError, 'Es muss Int[0]<=Int[1] sein'  
  
  #Feld-Dimension ermitteln
  s=Feld.shape
  if len(s)!=2:
    raise TypeError, 'Feld muss eine Matrix sein'
  nx=s[0]; ny=s[1]
  if nx%2!=1 or ny%2!=1:
    raise TypeError, 'Feld-Dimensionen muessen ungerade sein'
  
  #Int-Vektor in den k-Raum umrechnen
  #Arbeitskopie von Int erzeugen
  kInt=array(Int,float)
  if Art=='l':
    for i in xrange(2):
      if Int[i]==0:
        kInt[i]=inf
      elif isinf(Int[i]):
        kInt[i]=0
      else:
        kInt[i]=2*pi/Int[i]
    a=kInt[1]; kInt[1]=kInt[0]; kInt[0]=a
  #Das Quadrat von kInt
  kInt2=kInt**2.0
  
  #Filter notwendig
  if kInt[0]==0 and kInt[1]==inf:
    return Feld
  
  #dkxy ermitteln (fuer dkx ein mittleres dk)
  r_e=6371229.
  dkx=pi/(r_e*(nx-1)*deg2rad(dalpha)*cos((ny-1)*deg2rad(dalpha)/4.0))
  dky=pi/(r_e*(ny-1)*deg2rad(dalpha))
  
  #DCT
  f=dct2(Feld)
  
  #Filtern
  #Loop in C++-Code
  code = r"""
  double k2;
  for(int ix=0;ix<nx;ix++)
  {
    for(int iy=0;iy<ny;iy++)
    {
      k2=pow(ix*double(dkx),2)+pow(iy*double(dky),2);
      if((k2<kInt2(0))||(k2>kInt2(1))) f(ix,iy)=0.0;
    }
  }
  """
  err=weave.inline(code,['dkx','dky','nx','ny','kInt2','f'],type_converters=weave.converters.blitz,compiler='gcc',headers=["<math.h>"])
  
  #Loop in python
  #for ix in xrange(nx):
    #for iy in xrange(ny):
      ##Abstand zum Nullpunkt
      #k=sqrt((ix*dkx)**2+(iy*dky)**2)
      #if k<kInt[0] or k>kInt[1]:
	#f[ix,iy]=0.
  
  #IDCT
  return idct2(f)
  
def SpFilter2D(Feld,dx,Art,Int):
  
  #Wendet einen idealen spektralen Bandpassfilter auf das Feld mit dem aeuqidist. Gitterabstand dx an.
  #Art='k' oder 'l' gibt an, ob im Filterintervall Int k oder lambda gegeben ist. Die Einheit entspricht der von dx.
  #Int=[von bis] gibt das Intervall der Moden an, die beruecksichtigt werden sollen.
  
  #Eingabe ueberpruefen
  if Art!='l' and Art!='k':
    raise ValueError, 'Art muss k oder l sein'
  if len(Int)!=2:
    raise TypeError, 'Int muss 2dim Vektor sein'
  if Int[0]>Int[1]:
    raise ValueError, 'Es muss Int[0]<=Int[1] sein'  
  
  #Feld-Dimension ermitteln
  s=Feld.shape
  if len(s)!=2:
    raise TypeError, 'Feld muss eine Matrix sein'
  nx=s[0]; ny=s[1]
  if nx%2!=1 or ny%2!=1:
    raise TypeError, 'Feld-Dimensionen muessen ungerade sein'
  
  #Int-Vektor in den k-Raum umrechnen
  #Arbeitskopie von Int erzeugen
  kInt=array(Int,float)
  if Art=='l':
    for i in xrange(2):
      if Int[i]==0:
        kInt[i]=inf
      elif isinf(Int[i]):
        kInt[i]=0
      else:
        kInt[i]=2*pi/Int[i]
    a=kInt[1]; kInt[1]=kInt[0]; kInt[0]=a
  #dk ermitteln
  dkx=2*pi/(nx*dx)
  dky=2*pi/(ny*dx)
  
  #FFT und Feld um Nullpunkt zentrieren
  kFeld=fft.fftshift(fft.fft2(Feld))
  
  #Nullpunkt
  x0=(nx-1)//2
  y0=(ny-1)//2
  
  #Maske erzeugen
  Maske=zeros((nx,ny),float)
  #Iteration, erwuenschte Eintraege eins setzen
  #for ix in xrange(nx):
    #for iy in xrange(ny):
      ##Abstand zum Nullpunkt
      #k=sqrt(((ix-x0)*dkx)**2+((iy-y0)*dky)**2)
      #if k>=Int[0] and k<=Int[1]:
        #Maske[ix,iy]=1.
  
  code = r"""
  double k;
  for(int ix=0;ix<nx;ix++)
  {
    for(int iy=0;iy<ny;iy++)
    {   
      k=sqrt(pow((ix-x0)*dkx,2)+pow((iy-y0)*dky,2));
      if((k>=kInt(0))&&(k<=kInt(1))) Maske(ix,iy)=1;
    }
  }
  """
  err=weave.inline(code,['dkx','dky','x0','y0','nx','ny','kInt','Maske'],type_converters=weave.converters.blitz,compiler='gcc',headers=["<math.h>"])
  
  #Mittelwert erhalten oder filtern??
  #Hier: Mittelwert erhalten:
  Maske[x0,y0]=1.

  #Feld mit Maske multiplizieren
  kFeld*=Maske
  
  #Ruecktransformation
  kFeld=fft.ifft2(fft.ifftshift(kFeld))
  kFeld=kFeld.real
  
  return kFeld

def GaussLowPass2D(Feld,dx,Art,R):
  
  #Wendet einen Gauss-Lowpass-Filter mit Filterradius R auf das Feld mit dem aeuqidist. Gitterabstand dx an.
  #Art='k' oder 'l' gibt an, ob der Radius R als k oder lambda gegeben ist. Die Einheit entspricht der von dx.
  
  #Eingabe ueberpruefen
  if Art!='l' and Art!='k':
    raise ValueError, 'Art muss k oder l sein'
  if R<0:
    raise ValueError, 'R muss >=0 sein'  
  
  #Feld-Dimension ermitteln
  s=Feld.shape
  if len(s)!=2:
    raise TypeError, 'Feld muss eine Matrix sein'
  nx=s[0]; ny=s[1]
  if nx%2!=1 or ny%2!=1:
    raise TypeError, 'Feld-Dimensionen muessen ungerade sein'
  
  #R in den k-Raum umrechnen
  if Art=='l':
    R=2*pi/R
    
  #R quadrieren
  sqR=R*R
  
  #dk ermitteln
  dkx=2*pi/(nx*dx)
  dky=2*pi/(ny*dx)
  
  #FFT und Feld um Nullpunkt zentrieren
  kFeld=fft.fftshift(fft.fft2(Feld))
  
  #Nullpunkt
  x0=(nx-1)//2
  y0=(ny-1)//2
  
  #Filterfunkt erzeugen
  FF=zeros((nx,ny),float)
  
  code = r"""
  double sqk;
  for(int ix=0;ix<nx;ix++)
  {
    for(int iy=0;iy<ny;iy++)
    {   
      sqk=pow((ix-x0)*dkx,2)+pow((iy-y0)*dky,2);
      FF(ix,iy)=exp(-sqk/(2.0*sqR));
    }
  }
  """
  #err=weave.inline(code,['dkx','dky','x0','y0','nx','ny','sqR','FF'],type_converters=weave.converters.blitz,compiler='gcc',headers=["<math.h>"])
  
  for ix in xrange(nx):
    for iy in xrange(ny):
      sqk=((ix-x0)*dkx)**2+((iy-y0)*dky)**2
      FF[ix,iy]=exp(-sqk/(2.0*sqR))
  
  #Feld mit Filter multiplizieren
  kFeld*=FF
  
  #Ruecktransformation
  kFeld=fft.ifft2(fft.ifftshift(kFeld))
  kFeld=kFeld.real
  
  return kFeld
  

def dct(a,n=None,axis=-1):
  #Berechnet die DCT des Vektors a nach Typ 2 (Def. siehe wikipedia). Der Allgorithmus ist von
  #https://www.ee.washington.edu/techsite/papers/documents/UWEETR-2003-0026.pdf
  #n ist die Laenge des Vektors, der transformiert werden soll. a wird also entsprechend
  #gekuerzt oder mit Nullen gefuellt.
  #axis ist die Achse, entlang der die Transformation durchgefuert werden soll. axis=0 ist die Achse,
  #die durch Variieren des 0. Index entsteht.
  #Default fuer axis ist die letzte Achse
  dima=len(a.shape)
  if dima>2:
    raise TypeError, 'Max. 2-dim Felder werden unterstuetzt'
  if axis+1>dima:
    raise ValueError, 'axis zu gross'
  if axis==-1:
    axis=dima-1
  N=a.shape[axis]
  if n==None:
    n=N
  if dima==1:
    fz=zeros(2*n,float)
    if n>=N:
      fz[0:N]=a
    else:
      fz[0:n]=a[0:n]
    Y=fft.fft(fz)
    return real(Y[0:n]*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n)))
  if dima==2 and axis==0:
    fz=zeros((2*n,a.shape[1]),float)
    if n>=N:
      fz[0:N,:]=a
    else:
      fz[0:n,:]=a[0:n,:]
    Y=fft.fft(fz,axis=axis)
    return real((Y[0:n,:].transpose()*(exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n)))).transpose())
  if dima==2 and axis==1:
    fz=zeros((a.shape[0],2*n),float)
    if n>=N:
      fz[:,0:N]=a
    else:
      fz[:,0:n]=a[:,0:n]
    Y=fft.fft(fz,axis=axis)
    return real(Y[:,0:n]*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n)))
    
def idct(a,n=None,axis=-1):
  #Berechnet die IDCT des Vektors a nach Typ 2 (Def. siehe wikipedia). Der Allgorithmus ist von
  #https://www.ee.washington.edu/techsite/papers/documents/UWEETR-2003-0026.pdf
  #n ist die Laenge des Vektors, der transformiert werden soll. a wird also entsprechend
  #gekuerzt oder mit Nullen gefuellt.
  #axis ist die Achse, entlang der die Transformation durchgefuert werden soll. axis=0 ist die Achse,
  #die durch Variieren des 0. Index entsteht.
  #Default fuer axis ist die letzte Achse
  dima=len(a.shape)
  if dima>2:
    raise TypeError, 'Max. 2-dim Felder werden unterstuetzt'
  if axis+1>dima:
    raise ValueError, 'axis zu gross'
  if axis==-1:
    axis=dima-1
  N=a.shape[axis]
  if n==None:
    n=N
  if dima==1:
    fz=zeros(2*n,complex)
    if n>=N:
      fz[0:N]=a*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n))
    else:
      fz[0:n]=a[0:n]*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n))
    fz[0]*=1.0/n; fz[1:n]*=2.0/n
    Y=fft.fft(fz)
    return real(Y[0:n])
  if dima==2 and axis==0:
    fz=zeros((a.shape[1],2*n),complex)
    if n>=N:
      fz[:,0:N]=a.transpose()*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n))
    else:
      fz[:,0:n]=a[0:n,:].transpose()*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n))
    fz[:,0]*=1.0/n; fz[:,1:n]*=2.0/n
    Y=fft.fft(fz,axis=1)
    return real(Y[:,0:n].transpose())
  if dima==2 and axis==1:
    fz=zeros((a.shape[0],2*n),complex)
    if n>=N:
      fz[:,0:N]=a*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n))
    else:
      fz[:,0:n]=a[:,0:n]*exp(-complex(0,1)*array(xrange(n),float)*pi/(2.0*n))
    fz[:,0]*=1.0/n; fz[:,1:n]*=2.0/n
    Y=fft.fft(fz,axis=axis)
    return real(Y[:,0:n])

def testdct(f):
  #Berechnet die DCT des Vektors f nach Typ 2 (Def. siehe wikipedia) zu Testzwecken
  #Imlementiert die Definition und dient der Ueberpruefung
  N=f.size
  F=zeros(N,float)
  for k in xrange(N):
    for n in xrange(N):
      F[k]+=f[n]*cos(pi/N*(n+0.5)*k)
  return F
  
def testidct(f):
  #Berechnet die IDCT des Vektors f nach Typ 2 (Def. siehe wikipedia) zu Testzwecken
  #Imlementiert die Definition und dient der Ueberpruefung
  N=f.size
  F=ones(N,float)*f[0]*0.5
  for k in xrange(N):
    for n in xrange(1,N):
      F[k]+=f[n]*cos(pi/N*n*(k+0.5))
  return F*2.0/N
  
def dct2(a):
  #Berechnet die zweidimensionale dct
  dima=len(a.shape)
  if dima!=2:
    raise TypeError, 'Nur 2-dim Felder werden unterstuetzt'
  return dct(dct(a,axis=1),axis=0)
  
def idct2(a):
  #Berechnet die zweidimensionale dct
  dima=len(a.shape)
  if dima!=2:
    raise TypeError, 'Nur 2-dim Felder werden unterstuetzt'
  return idct(idct(a,axis=1),axis=0)
  
def detrend(a):
  #Entfernt den Trend aus a gemaess Errico: Spectra from LA-Grid
  dima=len(a.shape)
  if dima>2:
    raise TypeError, 'Max. 2-dim Felder werden unterstuetzt'
  if dima==1:
    print('spaeter!')
  if dima==2:
    nx=a.shape[0]; ny=a.shape[1]
    #tr speichert den Trend
    s=(a[nx-1,:]-a[0,:])/(nx-1)
    tr=outer(0.5*(2*seq(0,nx-1)-nx+1),s)
    da=a-tr
    t=(da[:,ny-1]-da[:,0])/(ny-1)
    tr+=outer(t,0.5*(2*seq(0,ny-1)-ny+1))
    da=a-tr
    return [da,tr]

def GeoWindp(phi,f,Lons,Lats):
  #Berechnet den geostrophischen Wind auf einer Druckflaeche
  #phi: Geopotenzial, f: Feld des Coriolisparameter, Lons, Lats: Koordinaten
  #Die Berechnung kann in geographischen oder rotierten Koordinaten durchgefuehrt werden. Die Ausgabe-Vektoren sind entsprechend definiert
  #Das Gitter muss regelmaessig sein. Die Randpunkte werden null gesetzt
  
  #Initialisierungen
  nlon=len(Lons);nlat=len(Lats); dx=Lons[1]-Lons[0]; dy=Lats[1]-Lats[0]; r_e=6371229.
  #Feld mit cos(Breite) erzeugen
  cosphi=ones((nlon,nlat),float)
  for ilat in xrange(nlat):
    cosphi[:,ilat]*=cos(deg2rad(Lats[ilat]))
  #Ausgabefelder anlegen
  ug=zeros((nlon,nlat),float); vg=zeros((nlon,nlat),float)
  #dphidlambda (nach Laenge ableiten)
  vg[1:nlon-1,:]=(phi[2:nlon,:]-phi[0:nlon-2,:])/(2.0*deg2rad(dx))
  #dphidphi (nach Breite ableiten)
  ug[:,1:nlat-1]=(phi[:,2:nlat]-phi[:,0:nlat-2])/(2.0*deg2rad(dy))
  #In Laengenableitungen umrechnen
  ug*=-1.0/(r_e*f)
  vg*=1.0/(r_e*f*cosphi)
  #Rueckgabe
  return [ug,vg]
  
def LeseTrajPos(File,nr):
  
  #Open file
  f = open(File, "r")
  buf=f.read(4)
  ntra=struct.unpack('>i',buf)[0]
  buf=f.read(4)
  ntrace=struct.unpack('>i',buf)[0]
  
  #result: x,y,z,tr1,tr2,...
  res=zeros(3+ntrace,float)
  
  lrec=4+(3+ntrace)*4
  #Search number
  i=1
  while(True):
    #Jump to position
    pos=8+(i-1)*lrec
    f.seek(pos,0)
    buf=f.read(4)
    if(buf==''):
      #eof
      x=nan; y=nan; z=nan
      break
    if(struct.unpack('>i',buf)[0]==nr):
      for j in xrange(3+ntrace):
        buf=f.read(4)
        res[j]=struct.unpack('>f',buf)[0]
      break
    i+=1
  
  return res
  
def LadeTraDef(Deffile):
  f = open(Deffile, "r")
  buf=f.read(4)
  ntra=struct.unpack('>i',buf)[0]
  buf=f.read(4)
  ntrace=struct.unpack('>i',buf)[0]
  DefM=zeros((ntra,5),float)
  lrec=20
  for i in xrange(ntra):
    pos=8+i*lrec
    f.seek(pos,0)
    buf=f.read(4)
    DefM[i,0]=struct.unpack('>i',buf)[0]
    buf=f.read(4)
    DefM[i,1]=struct.unpack('>f',buf)[0]
    buf=f.read(4)
    DefM[i,2]=struct.unpack('>f',buf)[0]
    buf=f.read(4)
    DefM[i,3]=struct.unpack('>f',buf)[0]
    buf=f.read(4)
    DefM[i,4]=struct.unpack('>i',buf)[0]
  f.close()
  return DefM
    
def FindeTrajNr(DefM,koord_in,tol=(0,0,0,0),koordsys='rot',npol=(0,90)):
  #Liefert die Trajektorien-Nummern mit den Koordinaten koord und der Toleranz tol zurueck
  koord=array(koord_in)
  if(koordsys=='geo'):
    koord[0:2]=geo2rot(koord[0:2],npol)
  matches=list()
  s=DefM.shape
  ntra=s[0]
  #Open Deffile
  for i in xrange(ntra):
    #Check
    #print(x,koord[0],y,koord[1],z,koord[2],t,koord[3])
    if((abs(DefM[i,1]-koord[0])<=tol[0]) and (abs(DefM[i,2]-koord[1])<=tol[1]) and (abs(DefM[i,3]-koord[2])<=tol[2]) and (abs(DefM[i,4]-koord[3])<=tol[3])):
      matches.append(int(DefM[i,0]))
  return matches

def geo2rot(k,npol):
  #rechnet geographische in rotierte Koordinaten um
  x=deg2rad(k[0]); y=deg2rad(k[1]); px=deg2rad(npol[0]); py=deg2rad(npol[1])
  rx=arctan(-cos(y)*sin(x-px)/(-cos(y)*sin(py)*cos(x-px)+sin(y)*cos(py)))
  ry=arcsin(sin(y)*sin(py)+cos(y)*cos(py)*cos(x-px))
  return (rad2deg(rx),rad2deg(ry))

def coarsegrain(Feld,n):
  #Coarse-Grained das Feld durch Mitteln ueber nxn Gitterpunkte.
  if(n==1): return Feld
  s=shape(Feld)
  dim=len(s)
  if(dim!=2):
    raise ValueError, 'Dim Feld muss 2 sein'
  nlon=s[0]; nlat=s[1]
  cnlon=nlon/n; cnlat=nlat/n
  
  cFeld=zeros((cnlon,cnlat),float)
  for ilon in xrange(cnlon):
    for ilat in xrange(cnlat):
      cFeld[ilon,ilat]=Feld[ilon*n:(ilon+1)*n,ilat*n:(ilat+1)*n].mean()
  return cFeld
  
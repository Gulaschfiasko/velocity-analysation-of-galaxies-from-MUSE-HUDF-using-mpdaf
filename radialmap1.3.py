import matplotlib.pyplot as plt  
import matplotlib # np und plt sind nur die kurzformen von 
import numpy as np                #das ist wie usepackage
import astropy.units as u
from mpdaf.obj import Cube        # das nimmt von einem package ne funktion damit wir nicht alles nehmen müssen oder so
from mpdaf.obj import Image, WCS 
from mpdaf.obj import gauss_image, moffat_image    # ich glaube das brauchen wir nichtmal
from scipy.interpolate import griddata #ich glaube das ist nur um gaussähnliche werte zu erstellen,
                                       # aber ich habe es vorsichtshalber trotzdem drin
from matplotlib.pyplot import figure, show     #haben wir das nicht schon mit dem obersten?
from matplotlib.ticker import MultipleLocator  #das ist auch für einen plot versuch
import lmfit
from lmfit.lineshapes import gaussian2d, lorentzian #das brauchen wir alles für 2D Gauss 
from mpdaf.sdetect.source import Source
from mpdaf.sdetect import Catalog
from scipy.optimize import curve_fit
import sys
import inspect

matplotlib.use("qt5agg")

np.seterr(divide='ignore', invalid='ignore')  #das ignoriert die fehlermeldungen 

ha = 6562.8 #wasserstoffline h_alpha in angström (656,28nm)
hb = 4861.327
o3 = 5006.9
c_light = 299792 #in km/s

glxy1 = Source.from_file('DR2_MOSAIC_000855.fits') #welches objekt wir uns ankucken
flux_filter = 9010
arz = 100 #ausreisserfilter
thicc = 0.5
#glxy1 = '124018038_objcube.fits'
#glxy1 = 'ADP.2019-09-24T07_25_24.757.fits'
#glxy1 = '102005037_objcube.fits'
#glxy1 = '101001006_objcube.fits'

tab = glxy1.tables['PL_Z']
z = tab[tab['FAMILY']==glxy1.REFZ]['Z'][0]

cube = glxy1.cubes['MUSE_CUBE']

#z = 0.24142        #  (kommt der aus dem umstellen der unteren gleichung? #0.24165 war drin, 0.24142 passt besser
                    # an sich sollten wir obs ja direkt messsen können)
                    # aus allgemeiner Datenbank ersichtlich (irgendwo online)
dl = 10.0            # delta lambda (fit fenster)
obs = ha * ( 1 + z ) #gleichung um rotverschiebung des observierten wertes herauszufinden
#wenn obs außerhalb des bereiches 
if obs > 9400.0:
    print("Linie außerhalb des Muse Spektrums:", obs)
    sys.exit(1)

print("Position der gewählten Emissionslinie bei:", obs , "Angström")

#cube = Cube(glxy1,ext=0) # erstellung eines cubes aus dem .fits file "glxy1"
                        #Was macht ext=0?(ansteuern der richtigen extension, wenn sie nicjt richtig benannt wurden)
cube.info()              # .info gibt alle zusätzlichen informationen die in der Datei gegeben sind aus
#cube sorted as z,y,x (also 0,1,2)
bild = cube.sum(axis=0) #wir definieren Bild als summe aller werte über die z achse,
                        # #also die helligkeiten aller bilder übereinander gelegt
"""
spektrum = cube.sum(axis=(1,2))
spektrum.plot()
plt.show()     #hier wird zusammen nach x und y summiert wodurch glaube ein spektrum entstehen müsste,
                # #man sieht also welche wellenlängen am meisten bedient werden
"""

#seg = bild.segment(background=500)
#source = seg[0]
#source.background()

#print (bild.data) #gibt riesige datenmatrix des bildes aus
print("********************")
print("2D-Gauss Fit Ergebnisse:") #automatische ausgabe der Werte eines Gauss Fits
gfit = bild.gauss_fit(plot=False) #der plot ist nicht umsetzbar, weil es anscheinend einen Fehler mit mpdaf gibt

a_elip, b_elip = gfit.fwhm
print(a_elip,b_elip)

gamma = gfit.rot
print("Rotation gegen Uhr von 12", gamma)

lz, ly, lx = cube.shape # definition der räumlichen ausdehnung

cntr = gfit.center
print("Die Mitte liegt bei ",cntr)

coordinates = bild.get_range()                                              #Hier berechnen wir die mitte des Gitters
print("Der Abschnitt geht über folgende Koordinaten ",coordinates)          #und rechnen sie in Bogensekunden um
print("Die y Ausdehnung ist ",coordinates[2]-coordinates[0])                #damit wir die genaue Position des Centers
print("Die x Ausdehnung ist ",coordinates[3]-coordinates[1])                #richtig auftragen können
x_int = (coordinates[3]-coordinates[1])/lx
print("x Intervall für ein Pixel", x_int, "in Grad bzw", x_int*3600, "In Bogensekunden" )
print("Die y Mitte des Gitters ist ", (coordinates[2]+coordinates[0])/2)
print("Die x Mitte des Gitters ist ", (coordinates[3]+coordinates[1])/2)

y_offset = ((coordinates[2]+coordinates[0])/2 - cntr[0])*3600 #dieser wert bei nutzung *-1 (kp warum)
x_offset = ((coordinates[3]+coordinates[1])/2 - cntr[1])*3600
print(y_offset)
print(x_offset)

#theta = gamma 

if a_elip < b_elip: # ja an sich geht auch b/a
    rela = b_elip/a_elip
    
else:
    rela = a_elip/b_elip
    
print( "Das Verhältnis beider Achsen ist 1 zu",rela)
# damit muss die Leuchtkraft erhöht werden
# abweichung durch dicke der scheibe und verteilung der sterne 
#(inklination muss keine verringerung der helligkeit bedeuten??)

# HIER MUSS ICH NOCHMAL MIT DER INKLINATION RICHTIG KUCKEN

"""
print("2D-Moffat Fit Ergebnisse")  #moffat wird als zuverlässiger für Teleskope am Erdboden betrachtet
mfit = bild.moffat_fit(plot=False) # Wendt hat gesagt Moffat ist unpraktisch
"""

lz, ly, lx = cube.shape # definition der räumlichen ausdehnung

bild.write('radmap.fits') # das bild was uns gleich gezeigt wird, wird im aktuellen verzeichnis gespeichert
bild.plot()
plt.show() #hier wird das vorher definierte Bild geplottet (also das bild der gesamten helligkeiten aller wellenlängen

#HIER BEGINNT DER TESTSTREIFEN

def f(x, y, A, x0, y0, sigma_x, sigma_y, theta):
    theta = np.radians(theta)
    sigx2 = sigma_x ** 2;
    sigy2 = sigma_y ** 2
    a = np.cos(theta) ** 2 / (2 * sigx2) + np.sin(theta) ** 2 / (2 * sigy2)
    b = np.sin(theta) ** 2 / (2 * sigx2) + np.cos(theta) ** 2 / (2 * sigy2)
    c = np.sin(2 * theta) / (4 * sigx2) - np.sin(2 * theta) / (4 * sigy2)

    expo = -a * (x - x0) ** 2 - b * (y - y0) ** 2 - 2 * c * (x - x0) * (y - y0)
    return A * np.exp(expo)

# N = 42                   # auflösung des meshgrids
x = np.linspace(-lx/10,lx/10,lx*10)  # x intervall (start, ende, pixel in x richtung)
y = np.linspace(-ly/10,ly/10,ly*10)  # in bogensekunden
# eigentlich war das dritte in der Klammer N und wir hatten die auflösung direkt eingestellt
# hier nurtezn wir jedoch auflösungsausgabe der cube.shape funktion

theta = gamma +90 # deg   #22  #rotation der galaxie
x0 = x_offset  #0.8             # xpos der galaxie
y0 = -y_offset  #0.5             # ypos der galxie
sigx = a_elip/2  #2             # a (breite der galaxie)
sigy = b_elip/2  #1             # b (höhe der galaxie)
A = 2          #2             # ich glaube das macht in unserem fall nichts mehr
Xg, Yg = np.meshgrid(x, y)  # meshgrid erzeugt quasi eine matrix mit den koordinaten x und y
Z = f(Xg, Yg, A, x0, y0, sigx, sigy, theta)

fig = figure()
ax = fig.add_subplot(1,1,1, aspect=1)
box = (x.min(), x.max(), y.min(), y.max())          # left, right, bottom, top
#statt Z hier überall jetzt bild.data, 28.03.
im = ax.imshow(bild.data, interpolation="none", origin='lower', extent=box, cmap='jet')

#print(Z) #gibt riesige datenmatrix von Z aus

ax.contour(Xg, Yg, Z , levels=[A/np.e**0.5], colors="k", extent=box) 
# erstellt eine ellipse um sigx und sigy (in der theorie)

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

ax.scatter(x0, y0)
ax.grid(color='w')
ax.xaxis.set_major_locator(MultipleLocator(1))

def rotate(vector, theta):
    th = np.radians(theta)
    v = np.array(vector)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    return R@v

xymaj = rotate((sigx,0), theta ) #das hier sorgt für den strich des a wertes der ellipse
xymaj2 = rotate((sigy,0), -1*theta )#hier versuche ich es auch für b zu replizieren
ax.plot( [x0,xymaj[0]+x0], [y0, xymaj[1]+y0] ), #ich glaube wir plotten hier eine art vektor vom definierten mittelpunkt zu einem anderen
ax.plot( [x0,xymaj2[1]+x0], [y0, xymaj2[0]+y0] ), #genau und hier ist nochmal der andere dazu

print(xymaj,xymaj2) #das gibt mir die einzelnen x und y werte der Vektoren aus

length_a = np.sqrt((xymaj[0])**2 + (xymaj[1])**2) #ausgabe betrag des Vektor xymaj 
print("Die große Halbachse (a) beträgt: ", length_a)

length_b = np.sqrt((xymaj2[0])**2 + (xymaj2[1])**2) #ausgabe betrag des Vektor xymaj 
print("Die kleine Halbachse (b) beträgt: ", length_b)




#HIER BEGINNT DIE ELLIPSENKORREKTUR
q_0 = 0.11  #dicke der galaxie bei betrachtung von edge on. Oft wird auch 0.075 genommen
height = q_0 * length_a #die (halbe) höhe einer galaxie  relativ zum Radius
b_approx = length_b
offset = 1
b_old = 3000 #zwei hohe werte, damit die Schleife nicht instant vorbei ist
b_new = 2000
 
while offset > 0.00001 :
    b_real = length_b - height * np.sin(np.arccos(b_approx / length_a))
    b_approx = b_real
    
    b_old = b_new
    b_new = b_real
    offset = 1 - b_new / b_old
    
    print(b_real)

#ELLIPSENKORREKTUR NACH WENDTS VORSCHLAG

q = length_b / length_a
inkl_alt = np.arccos(np.sqrt((q**2 - q_0**2)/(1 - q_0**2)))


print("Die angepasste kleine Halbachse beträgt:", b_real)
print("Das angepasste Verhältnis beträgt:", length_a/b_real)


inkl = np.arccos(length_b / length_a)
inkl_correct = np.arccos(b_real / length_a)

print("Die Ursprüngliche Inklination beträgt:", inkl, inkl * 180.0 / np.pi)
print("Die angepasste Inklination beträgt:", inkl_correct, inkl_correct * 180.0 / np.pi)
print("Die alternative Inklination beträgt:", inkl_alt, inkl_alt * 180.0 / np.pi)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

#flx_info = cube.flux
#print(flx_info)
flx = gfit.flux                                                    #Das hängt glaube auch von der Datei an
print("Gesamter Flux in 1e-20 erg / (Angstrom cm2 s):", flx)       #sich ab, aber ich weiß noch nicht ganz 
rflx = flx * length_a / length_b #realer Flux                      #wie ich das alleine auslese
print("Tatsächlicher Flux in 1e-20 erg / (Angstrom cm2 s):",rflx)  #Also wenn er zu uns gedreht ist
show()              # das zeigt den letzten Plot

allv=[] #eine leere Liste anlegen (allv steht nur für all values glaube ich)

print(lz,ly,lx) #ausgabe der räumlichen ausdehnung ³
print("Bitte Warten...")

for y in range(ly):#eine schleife 
    #print(y)       #ausgabe der Schleifennummer
    for x in range(lx):# noch eine schleife in der schleife
        spec = cube[:,y,x] #hier geben wir an was wir uns anschauen. Mit einem : kann man einen bereich kennzeichnen(0:4).
                            # Nur den : zu setzten bedeutet den gesamten Bereich
        #if y==20 and x>10 and x<30: #y ist für uns ungefähr die mitte, mit den beiden x geben wir den räumlichen bereich an in dem gesucht wird
         #   spec.plot()             #??? Warum Plotten wir hier 20 Bilder ???
          #  plt.show()
                
        try: #try ist ein befehl der nur sachen versucht, damit das program bei einem fehlschlag nicht direkt abstürzt
            testfit = spec.gauss_fit(lmin=obs-dl, lmax=obs+dl,plot=False) #wie liegt das pixel in dem bereich der Rotverschiebung
            center = testfit.lpeak  #Center einfach das center eines Punktes und nicht das Center der Galaxie?
            center_err=testfit.err_lpeak
            
            dv = (obs - center)/obs * c_light #berechnung der geschwindigkeit in km/s (rotverschiebung geeicht am galaxienkern)
        except:
            dv=np.nan # wenn der wert keine Zahl ist dann soll er rausgelassen werden
        #if center_err > 1.0: #or dv < -250 or dv > 250:
            #dv=np.nan # wenn Wert außerhalb des intervalls ist oder der Peak zu viele Pixel (1) vom Gausszentrum entfernt ist, dann soll der wert auch nicht berücksichtigt werden
        if bild[y,x]<flux_filter: #wenn der wert von der Lightmap unter diesem war, soll der dv wert nicht dargestellt werden
            dv=np.nan 
        else:
            allv.append(dv)    #alle werte die durchgekommen sind, kommen auf die liste
        bild[y,x]=dv           #auf das Pixel wird der dv wert gesetzt
        
        
v_max=max(allv)
v_min=min(allv)

if abs(v_max) > abs(v_min):
    v_top = abs(v_max)
else:
    v_top = abs(v_min)
        
print(np.average(allv), max(allv), min(allv)) #ausgabe des minimal und des maximalwertes (aus der liste) der rotation wegen ausreißern und den grenzen darunter immer mit vorschicht zu genießen
Imageanker = bild.plot(cmap=plt.cm.jet,vmin= -1*v_top,vmax=v_top) #jet ist die art und weise wie die farben dargestellt werden
cbar=plt.colorbar(Imageanker, cmap=plt.cm.jet) # 200 sind die grenzen von wo bis wo das farspektrum aus jet skaliert werden soll
cbar.set_label('radioal velocity [km/s]', rotation=270)
        
bild.write('colormap.fits')   
plt.title('Radial velocity map of Galaxy for Ha')
plt.savefig('plot_colormap.png')
plt.show()


allv.sort() #reverse=True
plt.plot(allv)
plt.show()

v_max = max(allv) / np.sin(inkl_alt)
v_min = min(allv) / np.sin(inkl_alt)

print("Die Tatsächlichen Rotationsgeschwindigkeiten sind", v_max, v_min)   
print("Im Durchschnitt ist die maximale Rotationsgeschwindigkeit der Galxie", (abs(v_max) + abs(v_min))/2)





bold = cube.sum(axis=0)

gfit = bold.gauss_fit(plot=False, unit_center=None) #das ist nur um nochmal sauber ans center zu kommen
cntr = gfit.center


print(gamma)
print(90-gamma)
steig = np.tan(-1*(90-gamma)/180*np.pi)
print("Steigung:", steig)
# if theta>45
# x_wert = coordinates[1] + x * x_int
#y_int = (coordinates[2]-coordinates[0])/ly #Bogensekunden abstand zweier pixel
#y_start = cntr[0] + (coordinates[1]-cntr[1])*steig  #in Winkel 
#y höhe der Achse x=0           
y_start = cntr[0] + (-cntr[1])*steig                   
#y_coord = (y_start - coordinates[0]) / y_int # y_start höhe in pixel
achse = []
#print("y_int=",y_int, y_int*3600)
print("y_start=",y_start,)
#print("y_coord=",y_coord)
for x in range(lx): #ist es nicht zu weit wenn er von 0 bis 42 zählt?
    for y in range(ly):
            spec = cube[:,y,x]
            y_wert = y_start  + x *  steig #funktion der achse (in pixel)
                                          
            
            etage = y #(coordinates[0] + y * y_int) #höhe der y ebene
            hoch = y_wert - etage    #abstand von der Funktion (nicht senkrecht)
            ds = hoch * np.cos(gamma/180*np.pi) #abweichung auf Funktion 
            dx = ds       #KEIN 90- DA HIN (auch wenn es dadurch besser aussieht uff)
            xs = x/np.sin((gamma)/180*np.pi) + dx   #tatsächlicher x wert schräg zur achse
            
            if bild[y,x] > v_min or bild[y,x] < v_max :
                #print("Der Bildpixel oder so", bild[y,x])
                dv = bild[y,x]         #dv wird wie bekannt genutzt
                achse.append([xs,dv])  #x und y werden zusammen auf eine Liste gemacht
 

achse.sort()            
#print(achse)    
achsenpixel = len(achse)  #Länge der Liste
pre_achspix = achsenpixel



print("x_int", x_int)

xs = [x[0] for x in achse] #beide werte werden wieder auf einzelne listen aufgeteilt
dv = [x[1] for x in achse] 
p=0

dv2=[] # neue liste wo alles rauf kommt, das nicht ausgefiltert wird

for p in range(achsenpixel): 
    q = p-1
    r = p-2
    xx=xs[p] 
    yy=dv[p]
    if p == 0 or p==1:
        dv2.append([xs[p],dv[p]])
    if p > 1:
        if abs (dv[p]) - abs(dv[r]) > -arz or abs(dv[p]) - abs(dv[q]) > -arz: #alles was sich im betrag um mehr als 100 in 
            dv2.append([xx,yy])                      #einem Schritt ändert kommt nicht mit auf die neue liste
                                                 #Funktioniert nicht bei Vorzeichenwechsel von 1. und 2. Wert
        
        
xs3 = [x[0] for x in dv2] #neue gefilterte liste wird für zweiten filter wieder aufgespalten
dv3 = [x[1] for x in dv2]

dv4=[]
        
for p in range(len(dv3)):
    q = p-1
    r = p-2
    if p==0 or p==1:
        dv4.append([xs3[p],dv3[p]])
    if p > 1:
        if dv3[r] - dv3[p] < arz or dv3[q] - dv3[p] < arz: #funktioniert nur wenn werte kleiner werden, aber auch bei
            dv4.append([xs3[p],dv3[p]])     #vorzeichenwechsel. Könnte aber zu overblocking führen

xs5 = [x[0] for x in dv4] 
dv5 = [x[1] for x in dv4] #/ np.sin(inkl_alt) #finale werte Wineklbereinigt

#print(xs5, dv5)

ausreisser = pre_achspix - len(dv5) #Werte auf liste vom anfang werden mit werte jetzt auf liste verglichen
print("Aussortierte Werte:", ausreisser)

for x in range(lx):         #wenn ein Pixel auf dem Weg zu dv5 aussortiert wurde
    for y in range(ly):     #wird er auch in bild[]aussortiert
        if bild[y,x] in dv5:
            lol=1
        else:
            bild[y,x] = np.nan

def tantan(x, k, a, m, b):
    return(np.arctan(x*k+a)*m+b)
    
xdata = np.asarray(xs5) #np.linspace(min(xs5), max(xs5), len(dv5))
ydata = np.asarray(dv5)

popt, pcov = curve_fit(tantan, xdata, ydata, p0=(0.3, -5, 100, -6)) #ungefähre parameter von 000855
popt                                                                # da sonst zu schnell mist raus kommt

#ax=plt.subplot(2, 1, 1)
print("Gefittete Werte:", popt)
dz = (popt[3]+np.average(dv5))/2 #geschwindigkeitsverschiebung um falschem z 
print(dz)                    # entgegen zu wirken
#tanvals = inspect.signature(tantan)
#print(tanvals)
#print("Abweichung des Nullpunktes:", b)
plt.plot(xdata, ydata)
plt.plot(xdata, tantan(xdata, *popt))
#plt.subplot(2, 1, 2,sharex=ax)
#plt.plot(xdata, tantan(xdata, *popt)-ydata, "x",c="tab:red")
plt.show()   

for x in range(lx):
    for y in range(ly):
        spec = cube[:,y,x] 
        bild[y,x] = bild[y,x] - dz
        
dv5 = dv5 - dz

v_max = max(dv5) #/ np.sin(inkl_alt)
v_min = min(dv5) #/ np.sin(inkl_alt)

if abs(v_max) > abs(v_min):
    v_top = abs(v_max)
else:
    v_top = abs(v_min)

print("Angepasst", np.average(dv5), max(dv5), min(dv5)) #ausgabe des minimal und des maximalwertes (aus der liste) der rotation wegen ausreißern und den grenzen darunter immer mit vorschicht zu genießen
Imageanker = bild.plot(cmap=plt.cm.jet,vmin=-1*v_top,vmax=v_top) #jet ist die art und weise wie die farben dargestellt werden
cbar=plt.colorbar(Imageanker, cmap=plt.cm.jet) # 200 sind die grenzen von wo bis wo das farspektrum aus jet skaliert werden soll
cbar.set_label('radioal velocity [km/s]', rotation=270)
        
bild.write('colormap.fits')   
plt.title('CORRECTED Z: Radial velocity map for Ha')
plt.savefig('plot_colormap.png')
plt.show()

plt.plot(xs) # überprüfung der äquidistanz (gerade=gut, treppen=schlecht)
plt.show()



for x in range(lx):#eine schleife 
    #print(y)       #ausgabe der Schleifennummer
    for y in range(ly):# noch eine schleife in der schleife
        spec = cube[:,y,x] #hier geben wir an was wir uns anschauen. Mit einem : kann man einen bereich kennzeichnen(0:4).
                            # Nur den : zu setzten bedeutet den gesamten Bereich
        y_wert = y_start  + x * steig
                     
        #thicc = 0.4    
                            
        if abs(y_wert - y) < thicc :
            dv=1000
            bild[y,x]=dv
        
        
print("Ursprünglich",np.average(allv), max(allv), min(allv)) #ausgabe des minimal und des maximalwertes (aus der liste) der rotation wegen ausreißern und den grenzen darunter immer mit vorschicht zu genießen
Imageanker = bild.plot(cmap=plt.cm.jet,vmin=-1*v_top,vmax=v_top) #jet ist die art und weise wie die farben dargestellt werden
cbar=plt.colorbar(Imageanker, cmap=plt.cm.jet) # 200 sind die grenzen von wo bis wo das farspektrum aus jet skaliert werden soll
cbar.set_label('radioal velocity [km/s]', rotation=270)
        
bild.write('colormap.fits')   
plt.title('WITH Main AXIS: Radial velocity map for Ha')
plt.savefig('plot_colormap.png')
plt.show()

print("Seitenverhältnis:", rela)
print("Inklination:", inkl_alt* 180.0 / np.pi)
print("v_top:", v_top)
print("echtes v_top:", v_top / np.sin(inkl_alt))
v_diff = max(dv5) - min(dv5)
print("v_diff:", v_diff)
print("Gesamter Flux in 1e-20 erg / (Angstrom cm2 s):", flx)     
print("Tatsächlicher Flux in 1e-20 erg / (Angstrom cm2 s):",rflx)
print("Rotverschiebung:", z)
tab = Catalog(glxy1.tables['HST_CAT'])
mag = tab['MAG_F775W'][0]
print("Magnitude:", mag)

v_new = np.average(dv5)
v_old = np.average(allv)
#print(v_new, v_old)

v_d = v_new - v_old
#print(v_d)
dz = z * v_d / c_light
#print(dz)
z_new = z-dz

print("Korrigiertes z:", z_new)

print(rela, inkl_alt* 180.0 / np.pi, v_top, v_top / np.sin(inkl_alt), v_diff, flx, rflx, z, mag, z_new)


#ALLE SPEKTREN AUF EIN GROßES AUFADDIEREN

v_curve = cube.sum(axis=0)
#gfit = v_curve.gauss_fit(plot=False)

v_curve.plot()
plt.show()

vdis = cube[:,int(ly/2),int(lx/2)]
vdis.plot()
plt.show()

lz, ly, lx = cube.shape

for y in range(ly):#eine schleife 
    #print(y)       #ausgabe der Schleifennummer
    for x in range(lx):# noch eine schleife in der schleife
        spec = cube[:,y,x] #hier geben wir an was wir uns anschauen. Mit einem : kann man einen bereich kennzeichnen(0:4).
                            # Nur den : zu setzten bedeutet den gesamten Bereich
        vdis = vdis + spec
        #if y==20 and x>10 and x<30: #y ist für uns ungefähr die mitte, mit den beiden x geben wir den räumlichen bereich an in dem gesucht wird
         #   spec.plot()             #??? Warum Plotten wir hier 20 Bilder ???
          #  plt.show()
                

vdis = vdis - cube[:,1,1]

vdis.plot()
plt.show()

#try: #try ist ein befehl der nur sachen versucht, damit das program bei einem fehlschlag nicht direkt abstürzt
#    testfit = vdis.gauss_fit(lmin=obs-dl*2, lmax=obs+dl*2,plot=True) #wie liegt das pixel in dem bereich der Rotverschiebung
#    center = testfit.lpeak  #Center einfach das center eines Punktes und nicht das Center der Galaxie?
#    center_err=testfit.err_lpeak
#    
#    dv = (obs - center)/obs * c_light #berechnung der geschwindigkeit in km/s (rotverschiebung geeicht am galaxienkern)
#
#    allv.append(dv)    #alle werte die durchgekommen sind, kommen auf die liste
#
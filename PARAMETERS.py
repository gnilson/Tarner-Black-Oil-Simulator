#Pressure at Standard Condictions, psia
psc = 14.7

#Initial Reservoir Pressure, psia
pi = 4000

#Bubble Point Pressure, psia
pb = 3000

#Pressure at Abandonment, psia
pab = 100

#Pressure Increments, psia
pstep = 100

#Oil Gravity, Deg. API
api = 25

#Gas Specific Graviy
SGgas = 0.7

#Reservoir Temperature, Deg. F
tres = 180

#Water Specific Gravity
SGwater = 1

#Relative Permeability to Gas @ Sor
krg_sor = 0.8

#Critical Gas Saturation
sgc = 0.05

#Irriducible Water Saturation
swi = 0.2

#Irriducible Oil Saturation
sor = 0.25

#Porosity
phi = 0.2

#Oil Formation Volume Factor at Bubble Point, RBBL/STB
bob = 1.277

#Solution Gas Oil Ratio at Bubble Point, SCF/STB
rsob = 500.3

#Formation Compressibility, 1/psi
cf = 5 * 10 ** -6

#Water Compessibility, 1/psi
cw = 3 * 10 ** -6

#Formation Thickness, ft
h = 50

#Drainage Area, ft^2
A = 40

#Initial Oil Formation Volume Factor, RBBL/STB
boi = 1.265

#Origional Oil In Place, STB
N = 7758 * A * h * phi * (1-swi) / boi

#Formation Permeability, mD
k = 100

#Pressure Drop Accross Wellbore, psi
deltapwf = 125

#Well Drainage Shape Factor
ca = 21.84

#Wellbore Radius, ft
rw = 0.5

#Skin Factor
sd = 0

mi = 0
vip = 0
pinj = 2500

#GOR Tolerance
tolerance = 0.0001

## read data sets
# China
# mortality data

## China male mortility data
data2 = read.table ('Chinamortalitymale.txt', header = F, sep = '')
data2pop = read.table ('Chinamortalitymalepop.txt', header = F, sep = '')
ages.mort = 0:90
years.mort = 1994:2010
China.mort.male = demogdata (data2, data2pop, ages.mort, years.mort, type = 'mortality', label = 'China', name = 'male')

## China female mortality data
data3 = read.table ('Chinamortalityfemale.txt', header = F, sep = '')
data3pop = read.table ('Chinamortalityfemalepop.txt', header = F, sep = '')
China.mort.female = demogdata (data3, data3pop, ages.mort, years.mort, type = 'mortality', label = 'China', name = 'female')

# China (presmooth)
China.mort.male.adjust = smooth.demogdata (China.mort.male, b = 65, k = 30)
China.mort.female.adjust = smooth.demogdata (China.mort.female, b = 65, k = 30)

China.lca.female = lca (China.mort.female.adjust, series = "female", adjust = "dt", max.age = 90,interpolate = TRUE)
ax.China.female = China.lca.female$ax
bx.China.female = China.lca.female$bx

China.lca.male = lca (China.mort.male.adjust, series = "male", adjust = "dt", max.age = 90, interpolate = TRUE)
ax.China.male = China.lca.male$ax
bx.China.male = China.lca.male$bx

kt.China.male = China.lca.male$kt
kt.China.female = China.lca.female$kt


# read multi-pop female mortality of 35 countries from Human Mortality Database
shortnames = c ("AUS","AUT","BLR","BGR","CAN","CHL","CZE","DNK","EST","FIN","FRATNP",
                "DEUTNP","HUN","ISL","IRL","ISR","ITA","JPN","LVA","LTU","LUX","NLD","NZL_NP",
                "NOR","POL","PRT","RUS","SVK","SVN","ESP","CHE","TWN","GBR_NP","USA","SWE")
names = c ("Australia","Austria","Belarus","Bulgaria","Canada","Chile","CzechRepublic",
           "Denmark","Estonia","Finland","France","Germany","Hungary","Iceland","Ireland","Israel",
           "Italy","Japan","Latvia","Lithuania","Luxembourg","Netherlands","NewZealand","Norway",
           "Poland","Portugal","Russia","Slovakia","Slovenia","Spain","Switzerland",
           "Taiwan","UnitedKingdom","USA","Sweden","China")

for (i in 1:35) {
  nam1 = paste (names [i] )
  assign (nam1, hmd.mx (shortnames[i], "fanglei@hu-berlin.de", "1440177160", names[i]))
  temp1 = hmd.mx (shortnames[i], "fanglei@hu-berlin.de", "1440177160", names[i])
  nam2 = paste (names[i], "lca.female", sep = ".")
  assign (nam2, lca (temp1, series = "female", adjust = "dt", interpolate = TRUE))
  temp2 = lca (temp1, series = "female", adjust = "dt", interpolate = TRUE)
  nam3 = paste ("ax", names[i], "female", sep = ".")
  assign (nam3, temp2$ax)
  nam4 = paste ("bx", names[i], "female", sep = ".")
  assign (nam4, temp2$bx)
  nam5 = paste ("kt", names[i], "female", sep = ".")
  assign (nam5, temp2$kt)
}
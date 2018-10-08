/*
The code contained in this file is a translation into Stata/Mata language of the original R code presented in the appendix of the following paper:

Rizzi, S., Halekoh, U., Thinggaard, M., Engholm, G., Christensen, N., Johannesen, T. B., & Lindahl-Jacobsen, R. (2018). **How to estimate mortality trends from grouped vital statistics**. _International journal of epidemiology._ [[PubMed]](https://www.ncbi.nlm.nih.gov/pubmed/30256946)
*/

version 14.2

cd "/Users/anddis/Dropbox/statistics/pclm2d/"
do "pclm2d.mata"

use "deaths.dta", clear
rename (Female Male Total) (d_=)
merge 1:1 Year Age using "exposures.dta", nogen
rename (Female Male Total) (e_=)
rename *, lower
keep if inrange(year, 1980, 2014)

// # exposures as offset
mata: E = st_data(., ("year", "age", "e_total"))

egen ageclass = cut(age), at(0(5)85 999)
foreach v of varlist ?_female ?_male ?_total {
	bysort year ageclass: egen sum_`v' = total(`v')
}
egen import_y = tag(year ageclass)

mata:
// # Deaths from cancer for both sexes combined
// # and put in vec form
y = st_data(., ("year", "ageclass", "sum_d_total"), "import_y")

// # Number of years
ny = length(uniqrows(y[,1]))

// # Number of bins
nb = J(1, ny, length(uniqrows(y[,2])))

// # Age grid for the underlying distribution for each calendar year
m = length(uniqrows(E[,2]))
x = 1..m

// # Define the grouping
// # e.g. 5 years age groups with 85+
ilo = uniqrows(y[,2])' :+ 1
ihi = ilo :+ 4
n = length(ihi)
ihi[n] = m
// # intervals lengths
leng = ihi :- ilo :+ 1

// # Construct C matrix
// # CA matrix in the age direction
CA = J(n, m, 0)
for (i=1; i<=n; i++) {
	CA[|i, ilo[i] \ i, ihi[i]|] = J(1, leng[i], 1)
}

// # CY matrix in the year direction
CY = diag(J(1, ny, 1))

// # C as kronecker product
Ci = CY # CA

// # exposures as offset
C = Ci * diag(E[, 3])

// # Construct B-spline basis 
// # for age
basisA = (1..m)'
xl = min(basisA)
xr = max(basisA)
xmax = xr + 0.01 * (xr - xl)
xmin = xl - 0.01 * (xr - xl)
BA = MortSmooth_bbase_R(basisA, xmin, xmax, floor(m/15), 3) 

// # for year
basisY = (1..ny)'
yl = min(basisY)
yr = max(basisY)
ymax = yr + 0.01 * (yr - yl)
ymin = yl - 0.01 * (yr - yl)
BY = MortSmooth_bbase_R(basisY, ymin, ymax, floor(ny/15), 3)
// # B as kronecker product
B = BY # BA

// # Second order penalties
DA = diff_R(I(cols(BA)), 2)
PA = I(cols(BY)) # (DA' * DA)
DY = diff_R(I(cols(BY)), 2)
PY = (DY' * DY) # I(cols(BA))

lambdaA_hat = 0.1
lambdaY_hat = 0.1
P = (lambdaA_hat :* PA) + (lambdaY_hat * PY)

// # Model estimation 
aTheta = aV = .
pclm2D(y[, 3], C, B, P, aTheta, aV)

log_gamma = B * aTheta
se_log_gamma = sqrt(rowsum(B :* (aV * B')'))

st_store(., st_addvar("double", ("log_gamma", "se_log_gamma")), (log_gamma, se_log_gamma))
end

gen lb_log_gamma = log_gamma - 1.96*se_log_gamma
gen ub_log_gamma = log_gamma + 1.96*se_log_gamma

gen obs_log_rate_sum = log(sum_d_total/sum_e_total)
gen obs_log_rate = log(d_total/e_total)

// figure 1
local brick ""
local a = 40
local col = 12
forv i = 1/7 {
	su log_gamma if age == `a' & year == 2014, meanonly
	local pos = r(mean)
	local brick `brick' (line log_gamma lb_log_gamma ub_log_gamma year if age == `a', lw(thin ..) ///
						lc(gs`col' ..) lp(l - -) text(`pos' 2018 "age `a'", col(gs`col'))) ///
						(scatter obs_log_rate year if age == `a', msiz(small) ms(Oh) mcol(gs`col'))
	local a = `a' + 10
	local col = `col' - 2
}
tw `brick', legend(off) xlabel(1980(5)2010 2014) xscale(range(1980 2020)) ///
	ytitle(log Mortality Rate) ylabel(, angle(horiz)) graphr(col(white) margin(1 1 1 1)) name(fig1, replace)
graph export "figure1.png", width(1000) replace
	
// figure 2
gen yearby = "{bf:" + strofreal(year) + "}" if inlist(year, 1980, 1990, 2000, 2010)
tw (scatter obs_log_rate age, msiz(small) ms(Oh) mcol(gs10)) ///
	(line log_gamma lb_log_gamma ub_log_gamma age, lc(black ..) lp(l - -) lw(thin ..)), ///
	 ytitle(log Mortality Rate) ylabel(, angle(horiz)) xtitle(Age) ///
	 xlabel(0(10)110, labsize(small)) ylabel(-10(2)0, labsize(small)) graphr(col(white)) ///
	 by(yearby, legend(off) graphr(col(white) margin(1 1 1 1)) note("") iy ix com) ///
	  subtitle(, bc(white)) name(fig2, replace)
drop yearby 
graph export "figure2.png", width(1000) replace

// figure 3
separate log_gamma if inrange(age, 50, 75), by(age)
separate obs_log_rate if inrange(age, 50, 75), by(age)

local lc gs14 gs12 gs8 gs6 gs0
local lp l - _ _._ dot 
local lw thin
tw (line obs_log_rate50-obs_log_rate59 year, lc(`lc' `lc') lp(`lp' `lp') lw(`lw' ..)) ///
(line obs_log_rate60-obs_log_rate69 year, lc(`lc' `lc') lp(`lp' `lp') lw(`lw' ..)) ///
(line obs_log_rate70-obs_log_rate75 year, lc(`lc' `lc') lp(`lp' `lp') lw(`lw' ..)), ///
	legend(off) xlabel(1980(5)2010 2014) ///
	ytitle(log Mortality Rate) ylabel(, angle(horiz)) graphr(col(white) margin(1 2 1 1)) name(fig3a, replace)

tw (line log_gamma50-log_gamma59 year, lc(`lc' `lc') lp(`lp' `lp') lw(`lw' ..)) ///
(line log_gamma60-log_gamma69 year, lc(`lc' `lc') lp(`lp' `lp') lw(`lw' ..)) ///
(line log_gamma70-log_gamma75 year, lc(`lc' `lc') lp(`lp' `lp') lw(`lw' ..)), ///
	legend(off) xlabel(1980(5)2010 2014) ///
	ytitle(log Mortality Rate) ylabel(, angle(horiz)) graphr(col(white) margin(1 2 1 1)) name(fig3b, replace)
	
graph combine  fig3a fig3b,  graphr(col(white) margin(1 2 1 1)) name(fig3, replace)
graph export "figure3.png", width(1000) replace
drop obs_log_rate50-obs_log_rate75 log_gamma50-log_gamma75

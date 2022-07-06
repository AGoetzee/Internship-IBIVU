# A single lazy script where I make all the free-energy surface plots for quick comparison.

library('metadynminer')
setwd('C:\\Users\\arthu\\Desktop\\simlab\\halfmeta')

hills_A <- read.hills('A\\HILLS_A')
hills_B <- read.hills('B\\HILLS_B')
hills_C <- read.hills('C\\HILLS_C')
hills_D <- read.hills('D\\HILLS_D')

A <- fes(hills_A)
B <- fes(hills_B)
C <- fes(hills_C)
D <- fes(hills_D)

plot(A,xlab='COM Separation (nm)',ylab='End to end distance (nm)',main='Scheme A')
plot(B,xlab='COM Separation (nm)',ylab='Radius of gyration (nm)',main='Scheme B')
plot(C,xlab='HisA:ThrI Distance (nm)',ylab='ThrA:HisI Distance (nm)',main='Scheme C')
plot(D,xlab='HisA:HisI Distance (nm)',ylab='ThrA:ThrI Distance (nm)',main='Scheme D')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\anglecv')
hills_angle_attach <- read.hills('attach\\HILLS')
angle_attach <- fes(hills_angle_attach)
plot(angle_attach, xlab='COM Separation (nm)',ylab='Angle (rad)', main='Scheme E')

hills_angle_sep <- read.hills('sep\\HILLS')
angle_sep <- fes(hills_angle_sep)
plot(angle_sep, xlab='COM Separation (nm)',ylab='Angle (rad)', main='Scheme E')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\compmeta\\80ns')
comp_hills <- read.hills('HILLS')
comp <- fes(comp_hills)
plot(comp,xlab='COM separation\nZ component (nm)',ylab='Relative FE (kJ/mol)',main='Distance components')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\1pdbmeta\\radgyr')
one_hills <- read.hills('HILLS')
one <- fes(one_hills)
plot(one,xlab='Radius of gyration (nm)',ylab='Relative FE (kJ/mol)',main='One chain')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\1pdbmeta\\e2e')
one_e2e_hills <- read.hills('HILLS')
one_e2e <- fes(one_e2e_hills)
plot(one_e2e,xlab='End-to-end distance (nm)',ylab='Relative FE (kJ/mol)',main='One chain')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\compmeta')
comp_C_hills <- read.hills('C\\HILLS')
comp_D_hills <- read.hills('D\\HILLS')
comp_D <- fes(comp_D_hills)
comp_C <- fes(comp_C_hills)

plot(comp_C,xlab='HisA:ThrI Z Distance (nm)',ylab='ThrA:HisI Z Distance (nm)',main='Scheme C')
plot(comp_D,xlab='HisA:HisI Z Distance (nm)',ylab='ThrA:ThrI Z Distance (nm)',main='Scheme D')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\boundaries\\500ns')
A_500ns_hills <- read.hills('A\\HILLS_A')
B_500ns_hills <- read.hills('B\\HILLS_B')
C_500ns_hills <- read.hills('C\\HILLS_C')

A_500ns <- fes(A_500ns_hills)
B_500ns <- fes(B_500ns_hills)
C_500ns <- fes(C_500ns_hills)

plot(A_500ns,xlab='COM Separation (nm)',ylab='End to end distance (nm)',main='Scheme A\n300 ns')
plot(B_500ns,xlab='COM Separation (nm)',ylab='Radius of gyration (nm)',main='Scheme B\n300 ns')
plot(C_500ns,xlab='HisA:ThrI Distance (nm)',ylab='ThrA:HisI Distance (nm)',main='Scheme C\n300ns')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\shortsep')
ss_A_hills <- read.hills('A\\HILLS_A')
ss_B_hills <- read.hills('B\\HILLS_B')
ss_C_hills <- read.hills('C\\HILLS_C')

ss_A <- fes(ss_A_hills)
ss_B <- fes(ss_B_hills)
ss_C <- fes(ss_C_hills)

plot(ss_A,xlab='COM Separation (nm)',ylab='End to end distance (nm)',main='Scheme A\n1.0 nm separation')
plot(ss_B,xlab='COM Separation (nm)',ylab='Radius of gyration (nm)',main='Scheme B\n1.0 nm separation')
plot(ss_C,xlab='HisA:ThrI Distance (nm)',ylab='ThrA:HisI Distance (nm)',main='Scheme C\n1.0 nm separation')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\combicv')
e2e_1_hills <- read.hills('1_e2e\\HILLS_1_e2e')
e2e_2_hills <- read.hills('2_e2e\\HILLS_2_e2e')
 
rg_2_hills <- read.hills('2_rg\\HILLS_2_rg')

e2e_1 <- fes(e2e_1_hills)
e2e_2 <- fes(e2e_2_hills)
rg_2 <- fes(rg_2_hills)

plot(e2e_1,ylab='End-to-end distance (nm)',xlab='Sum(distances)',main='Combination CV\nSum | E2E')
plot(e2e_2,ylab='End-to-end distance (nm)',xlab='Diff(distances)',main='Combination CV\nDifference | E2E')
plot(rg_2,ylab='Radius of gyration (nm)',xlab='Diff(distances)',main='Combination CV\nDifference | RadGyr')

setwd('C:\\Users\\arthu\\Desktop\\simlab\\rmsdcv\\target_meta')
rmsd_hills <- read.hills("HILLS")
rmsd <- fes(rmsd_hills,xlim=c(0,4))
rmsd_min <- fesminima(rmsd)
plot(rmsd)
plot(rmsd_min)

setwd('C:\\Users\\arthu\\Desktop\\simlab\\rmsdcv\\2dtarget')
rmsd_2d_hills <- read.hills('HILLS')
rmsd_2d <- fes(rmsd_2d_hills)
plot(rmsd_2d)
plot(fesminima(rmsd_2d))

setwd('C:\\Users\\arthu\\Desktop\\simlab\\rmsdcv\\')
run1_hills <- read.hills('run1\\HILLS')
run1 <- fes(run1_hills)
plot(run1,xlim=c(0,3),ylim=c(-5500,-4000))

interval.hills <- read.hills('interval\\HILLS_interval')
interval <- fes(interval.hills)
plot(interval)



setwd('C:\\Users\\arthu\\Desktop\\simlab\\rmsd-combi')
Crmsd.hills <- read.hills('HILLS_rmsd-c')
Crmsd <- fes(Crmsd.hills)
plot(Crmsd)
plot(fesminima(Crmsd))

cprmsd.hills <- read.hills('HILLS_rmsd-cp')
cprmsd <- fes(cprmsd.hills)
plot(cprmsd)
plot(fesminima(cprmsd))

e2e.rmsd.hills <- read.hills('HILLS_rmsd-e2e')
e2e.rmsd <- fes(e2e.rmsd.hills)
plot(e2e.rmsd)
plot(fesminima(e2e.rmsd))

setwd('C:\\Users\\arthu\\Desktop\\simlab\\restraints')
restraint.rmsd.hills <- read.hills('rmsd\\HILLS_rmsd')
restraint.rmsd <- fes(restraint.rmsd.hills)
plot(restraint.rmsd,xlim=c(0,4),ylim=c(-2100,-1900))


setwd('C:\\Users\\arthu\\Desktop\\simlab\\parabola')
parabola.hills <- read.hills('parabola\\HILLS_parabola')
par.restr.hills <- read.hills('par-restr\\HILLS_par-restr')
par.restr.rmsd.hills <- read.hills('par-restr-rmsd\\HILLS_par-restr-rmsd')
par.restr.rmsd2.hills <- read.hills('par-restr-rmsd-2\\HILLS_par-restr-rmsd-2')
par.params.hills <- read.hills('par-params\\HILLS_par-params')
par.uwall.hills <- read.hills('par-uwall\\HILLS_par-uwall')
par.newc.hills <- read.hills('par-new-c\\HILLS_par-new-c')
par.newc.nowall.hills <- read.hills('par-new-c-no-uwall\\HILLS_par-new-c-no-uwall')

parabola <- fes(parabola.hills)
par.restr <- fes(par.restr.hills)
par.restr.rmsd <- fes(par.restr.rmsd.hills)
par.restr.rmsd2 <- fes(par.restr.rmsd2.hills)
par.params <- fes(par.params.hills)
par.uwall <- fes(par.uwall.hills)
par.newc <- fes(par.newc.hills)
par.newc.nowall <- fes(par.newc.nowall.hills)

plot(parabola)
plot(par.restr)
plot(par.restr.rmsd)
plot(par.restr.rmsd2)
plot(par.params)
plot(par.uwall)
plot(par.newc)
plot(par.newc.nowall); plot(fesminima(par.newc.nowall)); linesonfes(neb(fesminima(par.newc.nowall),min1='A',min2='C'))

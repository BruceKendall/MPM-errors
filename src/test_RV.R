### test_RV.R
### Test out some ideas about reproductive value aggregation in unrolled models

library(mpmtools)

# A slightly simplified version of the lionfish lifetable
lf <- data.frame(stage = c("egg", "juvenile", "adult"),
                 survival = c(exp(-0.31*3) * exp(-0.35*27), exp(-0.165), exp(-0.052)),
                 maternity = c(0, 0, 194577), 
                 duration = c(1, 11, Inf))

A_unrolled_post <- make_stage4age_matrix(lf)
A_AAS_post <- make_stage4age_matrix(lf, approx_method = "AAS")

library(primer)
eigen_AAS_post <- DemoInfo(A_AAS_post)
eigen_unrolled_post <- DemoInfo(A_unrolled_post)

eigen_AAS_post$RV
eigen_unrolled_post$RV

eigen_AAS_post$SSD
sum(eigen_unrolled_post$SSD[2:12])


A_unrolled_pre <- make_stage4age_matrix(lf, model = "pre")
A_AAS_pre <- make_stage4age_matrix(lf, approx_method = "AAS", model = "pre")

eigen_AAS_pre <- DemoInfo(A_AAS_pre)
eigen_unrolled_pre <- DemoInfo(A_unrolled_pre)

eigen_AAS_pre$RV
eigen_unrolled_pre$RV

eigen_AAS_pre$SSD
sum(eigen_unrolled_pre$SSD[2:12])

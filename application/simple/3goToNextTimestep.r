source("../../efdm.r")

state0Area <- readRDS("./dat/state0Area.RData")
transition <- readRDS("./dat/transition.RData")

t1 <- efdmNextArea(state0Area, transition)
flow01Harv <- efdmNextFlow(state0Area, transition, "harvest")
flow01Resid <- efdmNextFlow(state0Area, transition, "residuals")

saveRDS(t1, file="./dat/t1.RData", compress="xz")
saveRDS(flow01Harv, file="./dat/harvest01.RData", compress="xz")
saveRDS(flow01Resid, file="./dat/residuals01.RData", compress="xz")


t2 <- efdmNextArea(t1, transition)
flow12Harv <- efdmNextFlow(t1, transition, "harvest")
flow12Resid <- efdmNextFlow(t1, transition, "residuals")
saveRDS(t2, file="./dat/t2.RData", compress="xz")
saveRDS(flow12Harv, file="./dat/harvest12.RData", compress="xz")
saveRDS(flow12Resid, file="./dat/residuals12.RData", compress="xz")

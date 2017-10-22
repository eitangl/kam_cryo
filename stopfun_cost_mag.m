function stopnow = stopfun_cost_mag(problem, x, info, last)

stopnow = info(last).cost < 1e-15;

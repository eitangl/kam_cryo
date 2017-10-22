function projections = gen_proj_ASPIRE(vol, K)

q = qrand(K);

projections = cryo_project(vol,q);
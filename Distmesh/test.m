fd=@(p) ddiff(drectangle(p,-1,1,-2,2),dcircle(p,0,0,0.5));
    fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
    [p,t]=distmesh2d(fd,fh,0.05,[-1,-2;1,2],[-1,-2;-1,2;1,-2;1,2]);
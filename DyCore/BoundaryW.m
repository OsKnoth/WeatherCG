function wCG=BoundaryW(v1CG,v2CG,CG,Param)
wCG=...
  -(Param.dXdxIC(:,:,:,1,3,1).*v1CG(:,:,:,1)...
  +Param.dXdxIC(:,:,:,1,3,2).*v2CG(:,:,:,1))...
  ./Param.dXdxIC(:,:,:,1,3,3);
end
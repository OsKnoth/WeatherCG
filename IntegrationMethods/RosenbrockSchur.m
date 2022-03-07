function V=RosenbrockSchur(V,dt,Fcn,Jac,CG,Param)
Vn=V;
ROS=Param.ROS;
nV1=size(V,1);
nV2=size(V,2);
nV3=size(V,3);
nJ=nV1*nV2*nV3;
nStage=ROS.nStage;
k=zeros([size(V) nStage]);
[JS]=Jac(V,CG,Param);
for iStage=1:nStage
  V=Vn;
  for jStage=1:iStage-1
    V=V+ROS.a(iStage,jStage)*k(:,:,:,jStage);
  end
  fV(:,:,:)=Fcn(V,CG,Param);
  for jStage=1:iStage-1
    fV=fV+(ROS.c(iStage,jStage)/(ROS.d*dt))*k(:,:,:,jStage);
  end
  k(:,:,:,iStage)=SchurSolve(fV,JS,dt*ROS.d);
end
V=Vn;
for iStage=1:nStage
  V=V+ROS.m(iStage)*k(:,:,:,iStage);
end
end

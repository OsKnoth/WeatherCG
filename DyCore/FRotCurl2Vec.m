function [RotCG]=FRotCurl2Vec(v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

vC1=reshape(...
  CG.DS*reshape(v2CG.*Param.dXdxIC(:,:,:,:,1,1)...
  -v1CG.*Param.dXdxIC(:,:,:,:,1,2)...
  ,OP,OP*NF*nz)...
  ,OP,OP,NF,nz)...
  -permute(...
  reshape(...
  CG.DS*reshape(...
  permute(...
  -v2CG.*Param.dXdxIC(:,:,:,:,2,1)...
  +v1CG.*Param.dXdxIC(:,:,:,:,2,2)...
  ,[2 1 3 4])...
  ,OP,OP*NF*nz)...
  ,OP,OP,NF,nz)...
  ,[2 1 3 4]);

D1cCG=reshape(...
  CG.DW*reshape(vC1,OP,OP*NF*nz)...
  ,OP,OP,NF,nz);
D2cCG=permute(reshape(...
  CG.DW...
  *reshape(...
  permute(...
  reshape(vC1,OP,OP,NF,nz)...
  ,[2 1 3 4])...
  ,OP,OP*NF*nz)...
  ,OP,OP,NF,nz)...
  ,[2 1 3 4]);
RotCG=zeros(OP,OP,NF,nz,2);
RotCG(:,:,:,:,2)=(-Param.dXdxIC(:,:,:,:,1,1).*D1cCG...
  -Param.dXdxIC(:,:,:,:,2,1).*D2cCG)...
  ./Param.JC;
RotCG(:,:,:,:,1)=(Param.dXdxIC(:,:,:,:,1,2).*D1cCG...
  +Param.dXdxIC(:,:,:,:,2,2).*D2cCG)...
  ./Param.JC;
end

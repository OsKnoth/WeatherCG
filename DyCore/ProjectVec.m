function [uS,vS]=ProjectVec(Fun,CG,Param)
OrdPoly=CG.OrdPoly;
OrdPolyZ=CG.OrdPolyZ;
nz=Param.Grid.nz;
uS=zeros(CG.NumG,nz,1);
vS=zeros(CG.NumG,nz,1);
uLoc3=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1);
vLoc3=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1);
for iz=1:nz
  for iF=1:Param.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        for k=1:OrdPolyZ+1
        x=Param.X(i,j,k,:,iF,iz);
        det=Param.J(i,j,k,iF,iz);
        [uLoc3(i,j,k),vLoc3(i,j,k)]=Fun(x,Param);
        uLoc3(i,j,k)=uLoc3(i,j,k)*det;
        vLoc3(i,j,k)=vLoc3(i,j,k)*det;
        end
      end
    end
    uLoc=0.5*(uLoc3(:,:,1)+uLoc3(:,:,2));
    vLoc=0.5*(vLoc3(:,:,1)+vLoc3(:,:,2));
    uS(CG.Faces(iF).Glob,iz)=uS(CG.Faces(iF).Glob,iz)+reshape(uLoc,(OrdPoly+1)*(OrdPoly+1),1);
    vS(CG.Faces(iF).Glob,iz)=vS(CG.Faces(iF).Glob,iz)+reshape(vLoc,(OrdPoly+1)*(OrdPoly+1),1);
  end
end
uS=uS./CG.M;
vS=vS./CG.M;
end

function [p]=Project(Fun,CG,Param)
OrdPoly=CG.OrdPoly;
OrdPolyZ=CG.OrdPolyZ;
nz=Param.Grid.nz;
p=zeros(CG.NumG,nz,1);
fLoc3=zeros(OrdPoly+1,OrdPoly+1);
for iz=1:nz
  for iF=1:Param.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        for k=1:OrdPolyZ+1
          x=Param.X(i,j,k,:,iF,iz);
          det=Param.J(i,j,k,iF,iz);
          fLoc3(i,j,k)=Fun(x,Param)*det;
        end
      end
    end
    fLoc=0.5*(fLoc3(:,:,1)+fLoc3(:,:,2));
    p(CG.Faces(iF).Glob,iz)=p(CG.Faces(iF).Glob,iz)+reshape(fLoc,(OrdPoly+1)*(OrdPoly+1),1);
  end
end
p=p./CG.M;
end

function ROS=RosenbrockMethod(Method)
switch Method
  case 'ROSEul'
    ROS.nStage=1;
    
  case 'ROS3Pw'
    ROS.nStage=3;
    ROS.beta0=7.88675134594812865529e-01;
    ROS.alpha(1,1)=2.0e0;
    ROS.alpha(2,1)=6.33974596215561403412e-01;
    ROS.alpha(3,1)=1.63397459621556140341e+00;
    ROS.alpha(3,2)=2.94228634059947813384e-01;
    ROS.alpha(3,3)=1.07179676972449078320e+00;
    
    ROS.beta(2,1)=-2.53589838486224561365e+00;
    ROS.beta(3,1)=-1.62740473580835520728e+00;
    ROS.beta(3,2)=-2.74519052838329002952e-01;
    ROS.alpha=ROS.alpha*ROS.beta0;
    ROS.beta=ROS.beta*ROS.beta0;
    ROS.d=ROS.beta0;
  case 'ROS2'
    ROS.nStage=2;
    % Matrix form aftrer Wolfbrandt
    %     d=.5;
    %     a(2,1)=2./3.;
    %     d(2,1)=.5;
    %     b(1)=.25;
    %     b(2)=.75
    alpha=zeros(ROS.nStage,ROS.nStage);
    gamma=zeros(ROS.nStage,ROS.nStage);
    ROS.d=0.5+sqrt(3)/6;
    d=ROS.d;
    alpha(2,1)=2/3;
    b(1,1)=0.25;
    b(1,2)=0.75;
    gamma(1,1)=d;
    gamma(2,1)=-d/b(1,2);
    gamma(2,2)=d;
    ROS.a=alpha/gamma;
    ROS.m=b/gamma;
    ROS.c=inv(diag(diag(gamma)))-inv(gamma);
  case 'ROSRK3'
    ROS.nStage=3;
    ROS.a=zeros(ROS.nStage,ROS.nStage);
    ROS.c=zeros(ROS.nStage,ROS.nStage);
    %beta0=1;
    %beta0=1.0/2.0
    beta0=(3.0+sqrt(3.0))/6.0;
    alpha(1,1)=1.0/(3.0*beta0);
    alpha(2,1)=(-1.0 + 12.0*beta0^2.0)/(18.0*beta0^2.0*(-1.0 + 4.0*beta0));
    alpha(2,2)=1.0/(2.0*beta0);
    alpha(3,1)=   (1.0 - 3.0*beta0*(7.0 + 16.0*beta0*(-2.0 + 3.0*beta0)))/...
      (36.0*beta0^3.0*(-1.0 + 4.0*beta0));
    alpha(3,2)=   (-1.0 + 12.0*beta0)/(4.0*beta0^2.0);
    alpha(3,3)=   1.0/beta0;
    beta(2,1)= (1.0 - 12.0*beta0^2.0)/(beta0^2.0*(-9.0 + 36.0*beta0));
    beta(3,1)= (-1.0 + 3.0*beta0*(7.0 + 16.0*beta0*(-2.0 + 3.0*beta0)))/...
      (36.*beta0^3.0*(-1.0 + 4.0*beta0));
    beta(3,2)= (0.250 - 3.0*beta0)/beta0^2.0;
    alpha=alpha*beta0;
    beta=beta*beta0;
    ROS.a(2:3,1:2)=alpha(1:2,1:2);
    ROS.m=alpha(3,:);
    ROS.c=beta;
    ROS.d=1;
end
end

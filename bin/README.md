# ./test_vanderpol

======================================
 >>   Backward Difference Formula    << 
 ======================================
  >> Number of variables    :    1
  >> Second order           : T
  >> Start time             :    0.000
  >> End time               :    2.000
  >> Step size              : 0.100E-02
  >> Number of steps        :       2001
  >> Physical System        : vanderpol 
  >> Max BDF Order          :    2
 Yet to deassociate the function
 fval         =   9.1728389816118288     
 Adjoint dfdx =   13.814951516273235     
 FD      dfdx =   13.814951749235592     
 Error        =   2.3296235696079748E-007
 Rel. Error   =   1.6863059762310623E-008

# ./test_pendulum

======================================
 >>   Backward Difference Formula    << 
 ======================================
  >> Number of variables    :    1
  >> Second order           : T
  >> Start time             :    0.000
  >> End time               :    3.140
  >> Step size              : 0.100E-02
  >> Number of steps        :       3141
  >> Physical System        : pendulum  
  >> Max BDF Order          :    2
  >> Opening file output/bdf.dat                  
 Yet to deassociate the function
 fval         =    1.5808405522020474     
 Adjoint dfdx =  -0.79016863004050120     
 FD      dfdx =  -0.79017049570140330     
 Error        =    1.8656609021006076E-006
 Rel. Error   =   -2.3610865151887681E-006

# ./test_adjoint

 ======================================
 >>   Backward Difference Formula    << 
 ======================================
  >> Number of variables    :    1
  >> Second order           : T
  >> Start time             :    0.000
  >> End time               :   10.000
  >> Step size              : 0.100E-02
  >> Number of steps        :      10001
  >> Physical System        : SMD       
  >> Max BDF Order          :    1
  >> Opening file output/bdf_forward_1.dat        
 Yet to deassociate the function
 fval         =   8.5906521075751598     
 Adjoint dfdx =   1.1419678851695612     
 FD      dfdx =   1.1419679044344861     
 Error        =   1.9264924988604548E-008
 Rel. Error   =   1.6869935585575612E-008

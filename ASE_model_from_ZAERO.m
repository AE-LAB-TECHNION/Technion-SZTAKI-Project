% Read ASE model at a certain airspeed that was generated as output of the ZAERO ASEOUT command
% It uses the sysrd.m function supplied by the ZAERO manual

clear all

filename='../ZAERO/ASE_PLANT_20'
% filename='../ZAERO/ASE_AE_20'
[V,rho,nx,nu,ny,ABCD,A,B,C,D,G,Bw,CG] = sysrd(filename); 
% nx number of states; nu - number of inputs; ny - number of outputs
% In the AE level: 4 modes (8 structural states) + 5 aero lags -> 13 states in Xae;
% In the PLANT level: 4 modes (8 structural states) + 5 aero lags -> 13 states in Xae; 4 CS with a 3rd order TF -> 12 actuator states; total nu=25 states

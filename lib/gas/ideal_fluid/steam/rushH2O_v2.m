function [r,u,s,h,cp,cv,a,g,Ts,ps]=rushH2O_v2(T,p,vapcalc);
%
% [r,u,s,h,cp,cv,a,Ts,ps]=rushH2O(T,p,vapcalc);
%
% Inputs
%     T       - temperature (K)
%     p       - pressure (Pa)
%     vapcalc - 0 for liquid
%               1 for superheated vapour
%               2 for metastable-vapour region
%
% Outputs
%     r       - density
%     u       - internal energy
%     s       - entropy
%     h       - enthalpy
%     cp      - constant pressure specific heat
%     cv      - constant volume specific heat
%     a       - speed of sound
%     g       - Gibbs free energy
%     Ts      - saturation temperature at p
%     ps      - saturation pressure at T
%
% Reference:
% Wagner et al., 2000, The IAPWS Industrial Formulation 1997 for 
% the Thermodynamic Properties of Water and Steam, 
% ASME J. Eng. Gas Turbines and Power, Vol 122, 150-182.
%
% David Buttsworth
% 10 Sept 2002 - version 1
% 22 July 2013 - version 2 - addition of metastable region data
%                          - something wrong with a for liquid state

% Gas constant
R=461.526;

% Calculate the saturation pressure and temperature
Tab=TableA11; ni=Tab(:,2);
% saturation p at given T
TT=T+ni(9)/(T-ni(10));
A=TT^2+ni(1)*TT+ni(2);
B=ni(3)*TT^2+ni(4)*TT+ni(5);
C=ni(6)*TT^2+ni(7)*TT+ni(8);
ps=1e6*(2*C/(-B+(B^2-4*A*C)^0.5))^4; 
% sat T at given p
PP=(p/1e6)^0.25;
E=PP^2+ni(3)*PP+ni(6);
F=ni(1)*PP^2+ni(4)*PP+ni(7);
G=ni(2)*PP^2+ni(5)*PP+ni(8);
D=2*G/(-F-(F^2-4*E*G)^0.5);
Ts=(ni(10)+D-((ni(10)+D)^2-4*(ni(9)+ni(10)*D))^0.5)/2; 

if vapcalc == 0,
    % Liquid state calculations
    Tab=Tab1; i=Tab(:,1); Ii=Tab(:,2); Ji=Tab(:,3); ni=Tab(:,4);
    tau=1386/T;
    pii=p/16.53e6;
    gamma=sum(ni.*(7.1-pii).^Ii.*(tau-1.222).^Ji);
    gamma_pi=sum(-ni.*Ii.*(7.1-pii).^(Ii-1).*(tau-1.222).^Ji);
    gamma_pipi=sum(ni.*Ii.*(Ii-1).*(7.1-pi).^(Ii-2).*(tau-1.222).^Ji);
    gamma_tau=sum(ni.*(7.1-pii).^Ii.*Ji.*(tau-1.222).^(Ji-1));
    gamma_tautau=sum(ni.*(7.1-pii).^Ii.*Ji.*(Ji-1).*(tau-1.222).^(Ji-2));
    gamma_pitau=sum(-ni.*Ii.*(7.1-pii).^(Ii-1).*Ji.*(tau-1.222).^(Ji-1));
    r=p/R/T/(pii*gamma_pi);
    u=R*T*(tau*gamma_tau-pii*gamma_pi);
    s=R*(tau*gamma_tau-gamma);
    h=R*T*tau*gamma_tau;
    cp=-tau^2*gamma_tautau*R;
    cv=R*(-tau^2*gamma_tautau+(gamma_pi-tau*gamma_pitau)^2/gamma_pipi);
    a=sqrt(R*T*gamma_pi^2/((gamma_pi-tau*gamma_pitau)^2/tau^2/gamma_tautau-gamma_pipi));
    g=gamma;
else % vapour state calculations
    if vapcalc == 1,
        Tab=TableA4; Ji0=Tab(:,2); ni0=Tab(:,3);
        Tab=TableA5; Ii=Tab(:,2); Ji=Tab(:,3); ni=Tab(:,4);
    elseif vapcalc == 2,
         Tab=TableA4_Eq23; Ji0=Tab(:,2); ni0=Tab(:,3);
         Tab=TableA7; Ii=Tab(:,2); Ji=Tab(:,3); ni=Tab(:,4);
    end
    tau=540/T;
    pii=p/1e6;
    gamma0=log(pii)+sum(ni0.*tau.^Ji0);
    gammar=sum(ni.*pii.^Ii.*(tau-0.5).^Ji);
    gamma=gamma0+gammar;
    gamma0_pi=1/pii;
    gammar_pi=sum(ni.*Ii.*pii.^(Ii-1).*(tau-0.5).^Ji);
    gamma0_tau=sum(ni0.*Ji0.*tau.^(Ji0-1));
    gammar_tau=sum(ni.*pii.^Ii.*Ji.*(tau-0.5).^(Ji-1));
    gamma0_pipi=-1/pii^2;
    gammar_pipi=sum(ni.*Ii.*(Ii-1).*pii.^(Ii-2).*(tau-0.5).^Ji);
    gamma0_tautau=sum(ni0.*Ji0.*(Ji0-1).*tau.^(Ji0-2));
    gammar_tautau=sum(ni.*pii.^Ii.*Ji.*(Ji-1).*(tau-0.5).^(Ji-2));
    gamma0_pitau=0;
    gammar_pitau=sum(ni.*Ii.*pii.^(Ii-1).*Ji.*(tau-0.5).^(Ji-1));
    r=p/R/T/(pii*(gamma0_pi+gammar_pi));
    u=R*T*(tau*(gamma0_tau+gammar_tau)-pii*(gamma0_pi+gammar_pi));
    s=R*(tau*(gamma0_tau+gammar_tau)-(gamma0+gammar));
    h=R*T*tau*(gamma0_tau+gammar_tau);
    cp=-R*tau^2*(gamma0_tautau+gammar_tautau);
    cv=-R*(tau^2*(gamma0_tautau+gammar_tautau)+(1+pii*gammar_pi-tau*pii*gammar_pitau)^2/(1-pii^2*gammar_pipi));
    a=sqrt(R*T*(1+2*pii*gammar_pi+pii^2*gammar_pi^2)/ ...
        ((1-pii^2*gammar_pipi)+ ...
        (1+pii*gammar_pi-tau*pii*gammar_pitau)^2/ ...
        (tau^2*(gamma0_tautau+gammar_tautau))));
    g=gamma*R*T;
end

function Tab=Tab1;
% Table 1: Liquid region ("region 1")
Table1=[
1 0 -2 0.146 329 712 131 67 +00
2 0 -1 -0.845 481 871 691 14 +00
3 0 0 -0.375 636 036 720 40 +01
4 0 1 0.338 551 691 683 85 +01
5 0 2 -0.957 919 633 878 72 +00
6 0 3 0.157 720 385 132 28 +00
7 0 4 -0.166 164 171 995 01 -01
8 0 5 0.812 146 299 835 68 -03
9 1 -9 0.283 190 801 238 04 -03
10 1 -7 -0.607 063 015 658 74 -03
11 1 -1 -0.189 900 682 184 19 -01
12 1 0 -0.325 297 487 705 05 -01
13 1 1 -0.218 417 171 754 14 -01
14 1 3 -0.528 383 579 699 30 -04
15 2 -3 -0.471 843 210 732 67 -03
16 2 0 -0.300 017 807 930 26 -03
17 2 1 0.476 613 939 069 87 -04
18 2 3 -0.441 418 453 308 46 -05
19 2 17 -0.726 949 962 975 94 -15
20 3 -4 -0.316 796 448 450 54 -04
21 3 0 -0.282 707 979 853 12 -05
22 3 6 -0.852 051 281 201 03 -09
23 4 -5 -0.224 252 819 080 00 -05
24 4 -2 -0.651 712 228 956 01 -06
25 4 10 -0.143 417 299 379 24 -12
26 5 -8 -0.405 169 968 601 17 -06
27 8 -11 -0.127 343 017 416 41 -08
28 8 -6 -0.174 248 712 306 34 -09
29 21 -29 -0.687 621 312 955 31 -18
30 23 -31 0.144 783 078 285 21 -19
31 29 -38 0.263 357 816 627 95 -22
32 30 -39 -0.119 476 226 400 71 -22
33 31 -40 0.182 280 945 814 04 -23
34 32 -41 -0.935 370 872 924 58 -25
];
i1=Table1(:,1);
Ii1=Table1(:,2);
Ji1=Table1(:,3);
ip=find(Table1(:,4)>0);
in=find(Table1(:,4)<0);
ni1(ip,1)=(Table1(ip,4)+Table1(ip,5)*1e-6+Table1(ip,6)*1e-9+Table1(ip,7)*1e-12+Table1(ip,8)*1e-14).*10.^Table1(ip,9);
ni1(in,1)=(Table1(in,4)-Table1(in,5)*1e-6-Table1(in,6)*1e-9-Table1(in,7)*1e-12-Table1(in,8)*1e-14).*10.^Table1(in,9);
Tab=[i1, Ii1, Ji1, ni1];

function Tab=Tab2;
% Table A1: boundary between regions 2 and 3
% Eq. (10) and (11)
Table2=[
1 0.348 051 856 289 69  +03
2 -0.116 718 598 799 75  +01
3 0.101 929 700 393 26  -02
4 0.572 544 598 627 46  +03
5 0.139 188 397 788 70  +02
];
i2=Table2(:,1);
ip=find(Table2(:,2)>0);
in=find(Table2(:,2)<0);
ni2(ip,1)=(Table2(ip,2)+Table2(ip,3)*1e-6+Table2(ip,4)*1e-9+Table2(ip,5)*1e-12+Table2(ip,6)*1e-14).*10.^Table2(ip,7);
ni2(in,1)=(Table2(in,2)-Table2(in,3)*1e-6-Table2(in,4)*1e-9-Table2(in,5)*1e-12-Table2(in,6)*1e-14).*10.^Table2(in,7);
Tab=[i2, ni2];

function Tab=TableA4;
% Table A4: Vapor region (ideal gas part) 
% Eq. (20) for use in Eq. (19)
Table=[
1 0 -0.969 276 865 002 17  +01
2 1 0.100 866 559 680 18  +02
3 -5 -0.560 879 112 830 20  -02
4 -4 0.714 527 380 814 55  -01
5 -3 -0.407 104 982 239 28  +00
6 -2 0.142 408 191 714 44  +01
7 -1 -0.438 395 113 194 50  +01
8 2 -0.284 086 324 607 72  +00
9 3 0.212 684 637 533 07  -01
];
i=Table(:,1);
Ji0=Table(:,2);
ip=find(Table(:,3)>0);
in=find(Table(:,3)<0);
ni0(ip,1)=(Table(ip,3)+Table(ip,4)*1e-6+Table(ip,5)*1e-9+Table(ip,6)*1e-12+Table(ip,7)*1e-14).*10.^Table(ip,8);
ni0(in,1)=(Table(in,3)-Table(in,4)*1e-6-Table(in,5)*1e-9-Table(in,6)*1e-12-Table(in,7)*1e-14).*10.^Table(in,8);
Tab=[i, Ji0, ni0];

function Tab=TableA4_Eq23;
% Table A4: Vapor region (ideal gas part) 
% Eq. (20) for use in Eq. (23) - metastable vapor region
Table=[
1 0 -0.969 372 683 930 49  +01      
2 1 0.100 872 759 700 06  +02
3 -5 -0.560 879 112 830 20  -02
4 -4 0.714 527 380 814 55  -01
5 -3 -0.407 104 982 239 28  +00
6 -2 0.142 408 191 714 44  +01
7 -1 -0.438 395 113 194 50  +01
8 2 -0.284 086 324 607 72  +00
9 3 0.212 684 637 533 07  -01
];
i=Table(:,1);
Ji0=Table(:,2);
ip=find(Table(:,3)>0);
in=find(Table(:,3)<0);
ni0(ip,1)=(Table(ip,3)+Table(ip,4)*1e-6+Table(ip,5)*1e-9+Table(ip,6)*1e-12+Table(ip,7)*1e-14).*10.^Table(ip,8);
ni0(in,1)=(Table(in,3)-Table(in,4)*1e-6-Table(in,5)*1e-9-Table(in,6)*1e-12-Table(in,7)*1e-14).*10.^Table(in,8);
Tab=[i, Ji0, ni0];

function Tab=TableA5;
% Table A5: Vapor region (residual part)
% Eq. (21) for use in Eq. (19)
Table=[
1 1 0 -0.177 317 424 732 13  -02
2 1 1 -0.178 348 622 923 58  -01
3 1 2 -0.459 960 136 963 65  -01
4 1 3 -0.575 812 590 834 32  -01
5 1 6 -0.503 252 787 279 30  -01
6 2 1 -0.330 326 416 702 03  -04
7 2 2 -0.189 489 875 163 15  -03
8 2 4 -0.393 927 772 433 55  -02
9 2 7 -0.437 972 956 505 73  -01
10 2 36 -0.266 745 479 140 87  -04
11 3 0 0.204 817 376 923 09  -07
12 3 1 0.438 706 672 844 35  -06
13 3 3 -0.322 776 772 385 70  -04
14 3 6 -0.150 339 245 421 48  -02
15 3 35 -0.406 682 535 626 49  -01
16 4 1 -0.788 473 095 593 67  -09
17 4 2 0.127 907 178 522 85  -07
18 4 3 0.482 253 727 185 07  -06
19 5 7 0.229 220 763 376 61  -05
20 6 3 -0.167 147 664 510 61  -10
21 6 16 -0.211 714 723 213 55  -02
22 6 35 -0.238 957 419 341 04  +02
23 7 0 -0.590 595 643 242 70  -17
24 7 11 -0.126 218 088 991 01  -05
25 7 25 -0.389 468 424 357 39  -01
26 8 8 0.112 562 113 604 59  -10
27 8 36 -0.823 113 408 979 98  +01
28 9 13 0.198 097 128 020 88  -07
29 10 4 0.104 069 652 101 74  -18
30 10 10 -0.102 347 470 959 29  -12
31 10 14 -0.100 181 793 795 11  -08
32 16 29 -0.808 829 086 469 85  -10
33 16 50 0.106 930 318 794 09  +00
34 18 57 -0.336 622 505 741 71  +00
35 20 20 0.891 858 453 554 21  -24
36 20 35 0.306 293 168 762 32  -12
37 20 48 -0.420 024 676 982 08  -05
38 21 21 -0.590 560 296 856 39  -25
39 22 53 0.378 269 476 134 57  -05
40 23 39 -0.127 686 089 346 81  -14
41 24 26 0.730 876 105 950 61  -28
42 24 40 0.554 147 153 507 78  -16
43 24 58 -0.943 697 072 412 10  -06
];
i=Table(:,1);
Ii=Table(:,2);
Ji=Table(:,3);
ip=find(Table(:,4)>0);
in=find(Table(:,4)<0);
ni(ip,1)=(Table(ip,4)+Table(ip,5)*1e-6+Table(ip,6)*1e-9+Table(ip,7)*1e-12+Table(ip,8)*1e-14).*10.^Table(ip,9);
ni(in,1)=(Table(in,4)-Table(in,5)*1e-6-Table(in,6)*1e-9-Table(in,7)*1e-12-Table(in,8)*1e-14).*10.^Table(in,9);
Tab=[i, Ii, Ji, ni];

function Tab=TableA7
% Table A7: Vapor region (residual part)
% Eq. (24) exponents for use in Eq. (23) - metastable vapor region
Table=[
1 1 0 -0.733 622 601 865 06 -02
2 1 2 -0.882 238 319 431 46 -01
3 1 5 -0.723 345 552 132 45 -01
4 1 11 -0.408 131 785 344 55 -02
5 2 1 0.200 978 033 802 07 -02
6 2 7 -0.530 459 218 986 42 -01
7 2 16 -0.761 904 090 869 70 -02
8 3 4 -0.634 980 376 573 13 -02
9 3 16 -0.860 430 930 285 88 -01
10 4 7 0.753 215 815 227 70 -02
11 4 10 -0.792 383 754 461 39 -02
12 5 9 -0.228 881 607 784 47 -03
13 5 10 -0.264 565 014 828 10 -02
];
i=Table(:,1);
Ii=Table(:,2);
Ji=Table(:,3);
ip=find(Table(:,4)>0);
in=find(Table(:,4)<0);
ni(ip,1)=(Table(ip,4)+Table(ip,5)*1e-6+Table(ip,6)*1e-9+Table(ip,7)*1e-12+Table(ip,8)*1e-14).*10.^Table(ip,9);
ni(in,1)=(Table(in,4)-Table(in,5)*1e-6-Table(in,6)*1e-9-Table(in,7)*1e-12-Table(in,8)*1e-14).*10.^Table(in,9);
Tab=[i, Ii, Ji, ni];



function Tab=TableA9;
% Table A9: 
% Eq. (25)
Table5=[
1 0 0 0.106 580 700 285 13  +01
2 0 0 -0.157 328 452 902 39  +02
3 0 1 0.209 443 969 743 07  +02
4 0 2 -0.768 677 078 787 16  +01
5 0 7 0.261 859 477 879 54  +01
6 0 10 -0.280 807 811 486 20  +01
7 0 12 0.120 533 696 965 17  +01
8 0 23 -0.845 668 128 125 02  -02
9 1 2 -0.126 543 154 777 14  +01
10 1 6 -0.115 244 078 066 81  +01
11 1 15 0.885 210 439 843 18  +00
12 1 17 -0.642 077 651 816 07  +00
13 2 0 0.384 934 601 866 71  +00
14 2 2 -0.852 147 088 242 06  +00
15 2 6 0.489 722 815 418 77  +01
16 2 7 -0.305 026 172 569 65  +01
17 2 22 0.394 205 368 791 54  -01
18 2 26 0.125 584 084 243 08  +00
19 3 0 -0.279 993 296 987 10  +00
20 3 2 0.138 997 995 694 60  +01
21 3 4 -0.201 899 150 235 70  +01
22 3 16 -0.821 476 371 739 63  -02
23 3 26 -0.475 960 357 349 23  +00
24 4 0 0.439 840 744 735 00  -01
25 4 2 -0.444 764 354 287 39  +00
26 4 4 0.905 720 707 197 33  +00
27 4 26 0.705 224 500 879 67  +00
28 5 1 0.107 705 126 263 32  +00
29 5 3 -0.329 136 232 589 54  +00
30 5 26 -0.508 710 620 411 58  +00
31 6 0 -0.221 754 008 730 96  -01
32 6 2 0.942 607 516 650 92  -01
33 6 26 0.164 362 784 479 61  +00
34 7 2 -0.135 033 722 413 48  -01
35 8 26 -0.148 343 453 524 72  -01
36 9 2 0.579 229 536 280 84  -03
37 9 26 0.323 089 047 037 11  -02
38 10 0 0.809 648 029 962 15  -04
39 10 1 -0.165 576 797 950 37  -03
40 11 26 -0.449 238 990 618 15  -04
];
i5=Table5(:,1);
Ii5=Table5(:,2);
Ji5=Table5(:,3);
ip=find(Table5(:,4)>0);
in=find(Table5(:,4)<0);
ni5(ip,1)=(Table5(ip,4)+Table5(ip,5)*1e-6+Table5(ip,6)*1e-9+Table5(ip,7)*1e-12+Table5(ip,8)*1e-14).*10.^Table5(ip,9);
ni5(in,1)=(Table5(in,4)-Table5(in,5)*1e-6-Table5(in,6)*1e-9-Table5(in,7)*1e-12-Table5(in,8)*1e-14).*10.^Table5(in,9);
Tab=[i5, Ii5, Ji5, ni5];

function Tab=TableA11;
% Table A11: Saturation pressure and temperature equations
% Coefficients of Eqs. (27), (28), and (55)
Table=[
1 0.116 705 214 527 67  +04
2 -0.724 213 167 032 06  +06
3 -0.170 738 469 400 92  +02
4 0.120 208 247 024 70  +05
5 -0.323 255 503 223 33  +07
6 0.149 151 086 135 30  +02
7 -0.482 326 573 615 91  +04
8 0.405 113 405 420 57  +06
9 -0.238 555 575 678 49  +00
10 0.650 175 348 447 98  +03
];
i=Table(:,1);
ip=find(Table(:,2)>0);
in=find(Table(:,2)<0);
ni(ip,1)=(Table(ip,2)+Table(ip,3)*1e-6+Table(ip,4)*1e-9+Table(ip,5)*1e-12+Table(ip,6)*1e-14).*10.^Table(ip,7);
ni(in,1)=(Table(in,2)-Table(in,3)*1e-6-Table(in,4)*1e-9-Table(in,5)*1e-12-Table(in,6)*1e-14).*10.^Table(in,7);
Tab=[i, ni];


Unit pj_startfinish;

interface

const numtypes=12;
      numbank=7;  {dont forget to change numbanks in this file}
      numtriggers=3;
      nummodules=3;
      nummulti=4;
      on=true;
      off=false;
      sideway=0;
      updown=1;
      totalnumber=84;
      maxnumofdifferenttrans=3;
      linedraw=-100;
      linedrawnoread=-101;
      repeatagain=1000;

type
     shortlines = array[1..26] of string[80];
     str30=string[130];   { increased from 30 (pj) }
     str15=string[15];
     vector=array[1..25] of str30;
     doubleint=array[1..2] of integer;
     intarray=array[1..25] of integer;
     int10array=array[1..10] of integer;
     str2010array=array[1..10] of string[20];
     str20=string[30];
     longvector=array[0..8191] of integer;
     bool84=array[1..84] of boolean;
     longintarray=array[0..999] of longint;
     int5array=array[1..5] of integer;
{     transspecarray=array[1..totalnumber,1..6] of string[11];}

var foreground,background,attribute:byte;
    clock:array[1..numbank] of byte;
    slope,coupling:array[1..numtriggers] of byte;
    d:integer;

    multiscnpage,graphmode,numbanks,xlow,ylow,jump,j,i,row,hotkey,numtrans,
    xcpos,ycpos,xlimitnumber,ylimitnumber,maxpix:integer;
    starttrans,finishtrans:array[1..maxnumofdifferenttrans] of integer;
    nummenu:array[0..30] of integer;
    transnum,channel:array[1..totalnumber] of integer;
    reclength:array[1..numbank] of integer;
    timebase:array[1..numbank,1..nummodules] of real;
    triggerunit,buffersize,tbase,pretrigsamp:array[1..3] of integer;
    start:array[1..numbank,1..nummodules,0..nummulti] of integer;
    dbaddressnum:array[1..5] of word;

    oldrow:doubleint;
    isize1:word;

    xtime:real;
    triglevel:array[1..numtriggers] of real;
    atset:array[1..numbank,1..nummodules] of real;
    qfluxgain,gain,sensitivity,location:array[1..totalnumber] of real;
    transducername,serial_number:array[1..totalnumber] of string[11];
    newxaxislimit,newyaxislimit:array[0..1] of str15;
    scale:array[1..10,1..14] of str15;
    qscale:array[1..6] of str15;
    ascale:array[1..7] of str15;
    fftscale:array[1..3] of str15;
    addscale:array[1..10] of str15;
    divisionscale,diffscale:array[1..11] of str15;
    fixscale:array[1..2] of str15;

    ch:char;

    imagepointer,textstart,textfinish:pointer;
    textmemory:^longintarray;
    dataarray:array[1..numbank,1..nummodules] of ^longvector;

    shortlineno : shortlines;
    filename,scramsurface,unitarray,zeroname,dbaddress:vector;
    currentfile:array[1..2] of string[15];
    destinationdisk:array[1..2] of str30;
    title:array[1..20] of string[50];
    monc_directory:string;

    graphinit,datapresent,modelconfigloaded,rundescriptionloaded,
    fftmade,drawn,quiet,errorsavingscreen,files_incomplete:boolean;
    dbchnl:array[1..numbank,1..nummodules,0..nummulti] of boolean;
    card:array[1..numbank] of boolean;

    { extra its for defrock }
    shotName, dataDir: string[130];
    withTimeColumn, physicalDataScales : boolean;


    procedure startup;
    procedure closedown;

implementation
{uses graph,messagefile,extras,crt,address;}

procedure startup;
var i,j,bank,module,multi,a,graphmode,graphdriver,result:integer;
    fname:text;
    dummy,nothing:str20;

begin
      assign(fname,'monc.val');
      {$I-}
      reset(fname);
      if ioresult=0 then
      begin
           for i:=1 to 25 do
           begin
                readln(fname,unitarray[i]);
                readln(fname,filename[i]);
                readln(fname,scramsurface[i]);
           end;
           unitarray[25]:=' ';
           {$V-}
           destinationdisk[1]:=filename[5];
           destinationdisk[2]:=filename[6];
           {$V+}
           for i:=1 to 2 do
           begin
                if length(destinationdisk[i])=1 then
                destinationdisk[i]:=destinationdisk[i]+':\';
                if copy(destinationdisk[i],length(destinationdisk[i]),1)<>'\'
                then destinationdisk[i]:=destinationdisk[i]+'\';
           end;
           for i:=1 to 10 do
           for j:=1 to 14 do readln(fname,scale[i,j]);
           for bank:=1 to numbank do
           for module:=1 to nummodules do
           for multi:=0 to nummulti do
           begin
                readln(fname,a);
                if a=1 then dbchnl[bank,module,multi]:=on
                else dbchnl[bank,module,multi]:=off;
           end;
           for i:=1 to 6 do readln(fname,qscale[i]);
           for i:=1 to 7 do readln(fname,ascale[i]);
           readln(fname,xcpos);
           readln(fname,ycpos);
           for i:=1 to 3 do readln(fname,fftscale[i]);
           for i:=1 to 10 do readln(fname,addscale[i]);
           for i:=1 to 11 do readln(fname,diffscale[i]);
           for i:=1 to 11 do readln(fname,divisionscale[i]);
           for i:=1 to 5 do readln(fname,zeroname[i]);
           for i:=1 to 5 do readln(fname,dbaddress[i]);
           for i:=1 to 2 do readln(fname,fixscale[i]);
           {$I+}
           close(fname);
      end
      {$I+}
      else
      begin
           writeln('Could not find monc.val in current directory.');
           for i:=1 to 25 do
           begin
                unitarray[i]:=' ';
                filename[i]:=' ';
                scramsurface[i]:=' ';
           end;
           for i:=1 to 10 do
           for j:=1 to 14 do scale[i,j]:='1';
           for bank:=1 to numbank do
           for module:=1 to 3 do
           begin
                start[bank,module,0]:=0;
                for multi:=0 to nummulti do dbchnl[bank,module,multi]:=on;
           end;
           for i:=1 to 6 do qscale[i]:=' ';
           for i:=1 to 7 do ascale[i]:=' ';
           for i:=1 to 3 do fftscale[i]:=' ';
           for i:=1 to 10 do addscale[i]:=' ';
           for i:=1 to 11 do diffscale[i]:='';
           for i:=1 to 11 do divisionscale[i]:='';
           for i:=1 to 5 do zeroname[i]:='';
           for i:=1 to 2 do fixscale[i]:='';
           dbaddress[1]:='320';
           dbaddress[2]:='310';
           dbaddress[3]:='311';
           dbaddress[4]:='e400';
           dbaddress[5]:='';
           for i:=1 to 4 do
           begin
                dummy:='$'+dbaddress[i];
                val(dummy,dbaddressnum[i],result);
           end;
           errorsavingscreen:=true;
      end;
      xlimitnumber:=0;
      ylimitnumber:=0;
      for i:=0 to 1 do
      begin
           newxaxislimit[i]:='0';
           newyaxislimit[i]:='0';
      end;
      for i:=1 to maxnumofdifferenttrans do
      begin
           starttrans[i]:=0;
           finishtrans[i]:=0;
      end;
end;

procedure closedown;
var i,j,bank,module,multi,a:integer;
    fname:text;
begin
    writeln('defrock finished');
end;

procedure initalizevalues;
var i,j:integer;
begin
     numbanks:=numbank;
     numtrans:=totalnumber;
     oldrow[1]:=0;
     nummenu[0]:=17;
     nummenu[1]:=3;
     nummenu[2]:=13;
     nummenu[3]:=5;
     nummenu[4]:=3;
     nummenu[5]:=9;
     nummenu[6]:=10;
     nummenu[7]:=8;
     nummenu[8]:=7;
     for i:=1 to numbanks do for j:=1 to nummodules do
     begin
          new(dataarray[i,j]);
          start[i,j,0]:=0;
     end;
     unitarray[25]:='';
     for i:=1 to 2 do currentfile[i]:=' No file ';
     for i:=1 to 25 do filename[i]:=' nothing ';
     quiet:=false;
     jump:=1;
     foreground:=1;
     background:=0;
     attribute:=16*background+foreground;
     for i:=1 to numbanks do reclength[i]:=2047;
     rundescriptionloaded:=false;
     modelconfigloaded:=false;
     datapresent:=false;
     graphinit:=false;
     drawn:=false;
     fftmade:=false;
     errorsavingscreen:=true;
     title[1]:='Raw data vs Time';
     title[2]:='Pressure vs Time';
     title[3]:='Thrust vs Time';
     title[4]:='Heat flux vs Time';
     title[5]:='Skin fric vs time';
     title[6]:='Massaged vs time';
     title[7]:='Pressure vs Distance';
     title[10]:='Magnitude vs Frequency';
     title[12]:='Norm press vs Time';
     title[13]:='Norm Thrust vs Time';
     title[14]:='Norm q dot vs time';
     title[15]:='Norm sk fr vs time';
     title[17]:='Norm Press vs Distance';
     title[18]:='Integrated Pressure vs Time';
     multiscnpage:=0;
end;

begin

     initalizevalues;
end.

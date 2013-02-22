unit pj_datadisk;
{$O+,F+}
interface
         procedure readmodel;  { moved from modeldis.pas }
         procedure readdata;
implementation
uses pj_startfinish;

var mstring:str20;
    startbank,endbank,startmodule,endmodule:integer;

procedure readmodel;
var fname:text;
    i,j:integer;
    mstring:str20;
    name:str30;
begin
     name:=filename[2]+'.mod';
     writeln('Reading Model File: ',name);
     if mstring<>'escape' then
     begin
          assign(fname,name);
          {$I-}
          reset(fname);
          if ioresult<>0 then
              begin
              writeln('Cannot find file: '+name);
              end
          else
              begin
              modelconfigloaded:=true;
              end;
          readln(fname,numtrans);
          for i:=1 to numtrans do
          begin
               readln(fname,transnum[i]);
               readln(fname,location[i]);
               readln(fname,sensitivity[i]);
               readln(fname,channel[i]);
               readln(fname,gain[i]);
          end;
          for i:=1 to 13 do readln(fname,scramsurface[i]);
          for i:=1 to numtrans do
          begin
               readln(fname,transducername[i]);
               if ioresult<>0 then transducername[i]:='';
          end;
          for i:=1 to numtrans do
          begin
               readln(fname,serial_number[i]);
               if ioresult<>0 then serial_number[i]:='';
          end;
          close(fname);
          {$I+}
     end;
end;


procedure readdata;
type str4=string[4];
var ext,mstring,directory:str20;
    name:str30;
    erroroccured:boolean;
    diskdrive:integer;


procedure readdatafile;
var bank,module:integer;
    messagestring:str20;
    fname:file;
procedure see_if_it_is_there;
begin
     str(bank*10+module,ext);
     if bank<10 then ext:='0'+ext;
     name:=filename[1]+'.'+ext;
     assign(fname,name);
     reset(fname);
     close(fname);
     if IOResult<>0 then
     begin
          writeln('cannot find file '+name);
          files_incomplete:=true;
     end;
end;

procedure briefcheck;
var found_lastone:boolean;
begin
     {$I-}
     files_incomplete:=false;
     found_lastone:=false;
     bank:=endbank+1;
     repeat
           module:=endmodule+1;
           bank:=bank-1;
           repeat
                 module:=module-1;
           until (dbchnl[bank,module,0]=on) or (module=startmodule);
           if dbchnl[bank,module,0]=on then found_lastone:=true;
     until found_lastone or (bank=startbank);
     if found_lastone then see_if_it_is_there;
end;


begin
     startbank:=1;
     endbank:=numbanks;
     startmodule:=1;
     endmodule:=nummodules;
     briefcheck;
     if files_incomplete=false then
     begin
          write('Read data files: ');
          for bank:=startbank to endbank do
          for module:=startmodule to endmodule do
          if dbchnl[bank,module,0]=on then
          begin
               str(bank*10+module,ext);
               if bank<10 then ext:='0'+ext;
               name:=filename[1]+'.'+ext;
               write(name+' ');
               {$I-}
               assign(fname,name);
               reset(fname,(reclength[bank]+1)*2);
               blockread(fname,dataarray[bank,module]^,1);{$I+}
               if ioresult<>0 then
                   writeln('Failed to read data from file: '+name)
               else
               begin
                   datapresent:=true;
                   currentfile[1]:=filename[1];
               end;
               close(fname);
          end;
     end;
end;

procedure read_header;
var realvar,state,a,b,c,d:real;
    bank,module,multi,i:integer;
    fname:file of real;
begin
     writeln('Reading ',filename[1],'.hed');
     name:=filename[1]+'.hed';
     assign(fname,name);
     {$I-}
     reset(fname);
     read(fname,realvar);
     numbanks:=round(realvar);
     if numbanks>numbank then numbanks:=numbank;
     for bank:=1 to numbanks do
        for module:=1 to nummodules do
           for multi:=0 to nummulti do
               begin
               read(fname,state);
               if state=0.
                   then dbchnl[bank,module,multi]:=on
                   else dbchnl[bank,module,multi]:=off;
               read(fname,realvar);
               start[bank,module,multi]:=round(realvar);
               end;
     for bank:=1 to numbanks do
          begin
          read(fname,realvar);
          reclength[bank]:=round(realvar);
          for module:=1 to nummodules do
               read(fname,timebase[bank,module]);
          end;

     for bank:=1 to numbanks do
     for module:=1 to nummodules do
         if dbchnl[bank,module,0]=on
             then read(fname,atset[bank,module])
             else atset[bank,module]:=5;
     for i:=1 to 3 do
         begin
         read(fname,a,b,c,d);
         triggerunit[i]:=round(a);
         buffersize[i]:=round(b);
         pretrigsamp[i]:=round(c);
         tbase[i]:=round(d);
         end;
     for i:=1 to numtriggers do
         begin
         read(fname,triglevel[i],a,b);
         slope[i]:=round(a);
         coupling[i]:=round(b);
         end;
     for i:=1 to numbanks do
         begin
         read(fname,c);
         clock[i]:=round(c);
         end;
     for i:=1 to totalnumber do qfluxgain[i]:=1;
     for i:=1 to totalnumber do read(fname,qfluxgain[i]);
     close(fname);
     {I+}
     if ioresult<>0 then writeln('Error reading hed file');
end;

procedure make_test_data;
var i,bank,module,multi:integer;
begin
          numbanks:=1;
          bank:=1;
          for module:= 1 to 3 do
          for multi:=0 to 4 do dbchnl[bank,module,multi]:=off;
          dbchnl[1,1,0]:=on;
          module:=1;
          multi:=0;
               start[bank,module,multi]:=1;
               reclength[bank]:=8196;
               timebase[bank,module]:=1;
     atset[bank,module]:=5;
     for i:=1 to 3 do
     begin
           triggerunit[i]:=1;
           buffersize[i]:=8192;
           pretrigsamp[i]:=10;
           tbase[i]:=1;
     end;
     for i:=1 to numtriggers do
     begin
          triglevel[i]:=0;
          slope[i]:=1;
          coupling[i]:=1;
     end;
     for i:=1 to numbanks do
     begin
          clock[i]:=1;
     end;
     for i:=1 to totalnumber do qfluxgain[i]:=1;
     for i:=1 to 200 do dataarray[bank,module]^[i]:=0;
     for i:=201 to 8192 do dataarray[bank,module]^[i]:=round(i/2);
     datapresent:=true;
end;

begin               {readdata}
     diskdrive:=1;
     if filename[1]='test' then make_test_data
     else
     begin
          writeln('Begin readdata...');
          writeln('call read_header');
          read_header;
          writeln('call readdatafile');
          readdatafile;
     end;
     writeln('End of readdata.');
end;
end.
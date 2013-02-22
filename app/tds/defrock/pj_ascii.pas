unit pj_asciiconvert;
{$O+,F+}
interface
         procedure writetodisc;
implementation
uses pj_startfinish;

var menuarray,scalearray:^vector;
    numofinputs,numsigfig:integer;
    datatype:intarray;
    lngth:doubleint;
    textfile:text;
    timewanted:boolean;

procedure writetodisc;
var bank,module,multi,startnum,finishnum,numofdpts,numofavpts,totnumavpts,
    minn,maxn,i,j,result,stpt:integer;
    maxtime,mintime,averagetime,timescale,yscale,y,t:real;
    error,integrated:boolean;
    ext,mstring,directory:str20;
    timeInterval, timeStart : real;
    unitIndex : integer;


procedure splitchannel(i:integer);
begin
     bank:=trunc(channel[i]/100);
     module:=abs(trunc((channel[i]-bank*100)/10));
     multi:=abs(trunc(channel[i]-bank*100-module*10));
     bank:=(bank-1) mod 7 +1;
     module:=(module-1) mod 3+ 1;
     multi:=multi mod 5;
end;

function zeropt(bank,module,multi:integer):real;
var i:integer;
    zero:real;
begin
     zero:=0;
     for i:=1 to 20 do zero:=zero+dataarray[bank,module]^[i+start[bank,module,multi]];
     zeropt:=zero/20;
end;

procedure gettimescale;
var firstch:string[1];
begin
     firstch:=copy(unitarray[7],1,1);
     if (firstch='u') or (firstch='U') then timescale:=1
     else
     if (firstch='s') or (firstch='S') then timescale:=1e6
     else timescale:=1e3;
end;

procedure getyscale(chnlnum,bank,module:integer);
begin
    if physicalDataScales then
        begin { Allan's complete scaling }
        yscale:=atset[bank,module]/2048/sensitivity[chnlnum]/gain[chnlnum]/qfluxgain[chnlnum];
        if integrated then
            if ((transnum[chnlnum]=3) or (transnum[chnlnum]=4)) and
            ((copy(unitarray[3],1,1)='M') or (copy(unitarray[3],1,1)='m')) then
            yscale:=yscale*1e-6;
        end
    else
        begin { Scale only to recorded voltage. }
        yscale:=atset[bank,module]/2048/gain[chnlnum]/qfluxgain[chnlnum];
        end
end;

procedure getextremes(bank,module:integer);
var a:integer;
procedure getlengthofrecord;
var i,number:integer;
begin
     number:=0;
     for i:=1 to nummulti do
     if dbchnl[bank,module,i]=on then number:=number+1;
     if number=0 then
          numofdpts:=reclength[bank]
     else numofdpts:=trunc(reclength[bank]/number);
end;
begin      {getextremes}
     getlengthofrecord;
     numofavpts:=trunc(averagetime*timescale/timebase[bank,module]/2);
     if numofavpts>numofdpts/2 then numofavpts:=trunc(numofdpts/2);
     totnumavpts:=2*numofavpts+1;

     {take the lot }
     minn := numofavpts;
     maxn := numofdpts;
     { minn:=trunc(mintime*timescale/timebase[bank,module]); }
     { maxn:=trunc(maxtime*timescale/timebase[bank,module]); }

     if minn<numofavpts then minn:=numofavpts;
     if minn>numofdpts then minn:=numofdpts;
     if maxn<numofavpts then maxn:=numofavpts;
     if maxn>numofdpts then maxn:=numofdpts;
     if minn>maxn then
     begin
          a:=minn;
          minn:=maxn;
          maxn:=a;
     end;
end;


procedure writefirstpt(var y:real;locationnum:integer);
var j:integer;
begin
     for j:=minn-numofavpts to minn+numofavpts do
         y:=y+dataarray[bank,module]^[stpt+j];
     if physicalDataScales then
         y:=yscale*(y/totnumavpts-zeropt(bank,module,multi))
     else
         y:=yscale*(y/totnumavpts);
     t:=1000.0*minn*timebase[bank,module]/timescale;  {microsec}
     if timewanted then
         writeln(textfile,t:12:-5,' ',y:12:-5)
     else
         writeln(textfile,y:12:-5);
end;

begin  {writetodisc}
     if datapresent then
     if modelconfigloaded then
     begin
          error:=false;
          startnum:=1;  { do all channels }
          finishnum:=numtrans;
          if not error then
          begin
               gettimescale;
               averagetime := 0.0;
               timewanted := withTimeColumn;
               integrated:=false;
               directory:='';
               write('Writing ASCII files: ');
               for i:=startnum to finishnum do
               begin
                    str(channel[i],ext);
                    if mstring<>'escape' then
                    begin
                         error:=false;
                         assign(textfile,filename[4]+'.'+ext);
                         {$I-}
                         rewrite(textfile);
                         if ioresult<>0 then error:=true
                         else
                         begin
                              write(filename[4],'.',ext,' ');
                              splitchannel(i);
                              getyscale(i,bank,module);
                              getextremes(bank,module);
                              stpt:=start[bank,module,multi];

                              {start of textfile header }
                              writeln(textfile, '# dataSource MONC v4.8, extracted by defrock');
                              writeln(textfile, '# shotName ', shotName);
                              writeln(textfile, '# channelId ', channel[i]);
                              if withTimeColumn then
                                  writeln(textfile, '# withTimeColumn yes')
                              else
                                  writeln(textfile, '# withTimeColumn no');
                              writeln(textfile, '# dataPoints ', maxn - minn + 1);
                              case transnum[i] of
                                  1,2,7 : unitIndex := 1;
                                  3,4   : unitIndex := 3;
                                  5,6   : unitIndex := 5;
                                  9,10  : unitIndex := 6;
                                  11    : unitIndex := 8;
                                  else unitIndex := 25
                              end;
                              if physicalDataScales then
                                  begin
                                  writeln(textfile, '# dataType scaled');
                                  writeln(textfile, '# dataUnits ', unitarray[unitIndex]);
                                  end
                              else
                                  begin
                                  writeln(textfile, '# dataType raw');
                                  writeln(textfile, '# dataUnits volts')
                                  end;
                              timeInterval := 1000.0*timebase[bank,module]/timescale; {microsec}
                              timeStart := timeInterval * minn;
                              writeln(textfile, '# timeStart ', timeStart);
                              writeln(textfile, '# timeAverageWindow ', averagetime);
                              writeln(textfile, '# timeInterval ', timeInterval);
                              writeln(textfile, '# timeUnits microseconds');
                              writeln(textfile, '# transducerSensitivity ', sensitivity[i]);
                              writeln(textfile, '# transducerSensitivityUnits volts/', unitarray[unitIndex]);
                              writeln(textfile, '# transducerName ', transducername[i]);
                              writeln(textfile, '# transducerLocation ', location[i]);
                              writeln(textfile, '# transducerSerialNumber ', serial_number[i]);
                              writeln(textfile, '# gain ', gain[i]);
                              writeln(textfile, '# qfluxgain ', qfluxgain[i]);
                              writeln(textfile, '');
                              { end of textfile header }

                              y:=0;
                              writefirstpt(y,i);
                              if ioresult<>0 then error:=true else
                              begin
                                   j:=minn;
                                   if minn<maxn then
                                   repeat
                                         j:=j+1;
                                         t:=1000.0*j*timebase[bank,module]/timescale;{microsec}
                                         y:=y+(dataarray[bank,module]^[j+stpt+numofavpts]-
                                         dataarray[bank,module]^[stpt+j-1-numofavpts])*yscale/totnumavpts;
                                         if timewanted then
                                             writeln(textfile,t:12:-5,' ',y:12:-5)
                                         else
                                             writeln(textfile,y:12:-5);
                                         if ioresult<>0 then error:=true;
                                   until (error=true) or (j=maxn);
                                   close(textfile);
                                   {$I+}
                              end;
                         end;
                         if error=true then writeln('Disc is full or stuffed.');
                    end;
               end;
          end else writeln('Channel number '+ext+' not recorded.');
     end else writeln('Model configuration not present')
     else writeln('Data is not present');

     writeln('');
     writeln('Write list of channels and transducer names');
     assign(textfile,filename[4]+'.lst');
     rewrite(textfile);
     for i:=startnum to finishnum do
         begin
         writeln(textfile, channel[i], ' ', transducername[i]);
         end;
     close(textfile);

     writeln('end of writetodisc');
end;



end.


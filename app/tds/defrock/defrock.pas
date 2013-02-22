program defrock;

(* Convert recently generated MONC binary files to text.
 * This program was adapted from the MONC source files,
 * specifically startfin.pas, datadisk.pas, modeldis.pas and asciicon.pas.
 *
 * PJ, 27-dec-01
 *     13-Jan-02: Adapted to Kylix on Linux
 *)


{$F+}
{$REALCOMPATIBILITY ON}

uses pj_startfinish, pj_datadisk, pj_asciiconvert;

var homeDir: string[130];

begin
     writeln('Begin defrock...');

     { Defaults }
     withTimeColumn := false;
     physicalDataScales := true;

     GetDir(0, homeDir);  { remember where we started }

     { Normally, we shall expect at least 2 command-line parameters. }
     if ParamCount < 2 then
     begin
         writeln('');
         writeln('Usage: defrock dataDir shotName [notimes|withtimes [volts]]');
         writeln('');
         writeln('Required command-line parameters:');
         writeln('dataDir is the path from which one can see the directory shotName');
         writeln('which contains the individual binary data files.');
         writeln('This directory should also contain the rundesc and descript directories.');
         writeln('The ASCII (text) files will be written to the same directory as ');
         writeln('the original binary data files but will have names shotNameA.nnn ');
         writeln('where nnn is the channel number.');
         writeln('By default, defrock will write only a sample value column ');
         writeln('and the sample values will be in physical units.');
         writeln('Optional command-line parameters:');
         writeln('The parameter "withtimes" will cause the sample times to be included while');
         writeln('"notimes" will suppress the time column. (default="notimes")');
         writeln('The parameter "volts" will cause the sample values to be left as voltages.');
         writeln('');
         writeln('Example 1 : minimal data file with data already scaled (use / on Unix)');
         writeln('defrock data\T4\ 7319');
         writeln('Example 2 : files ready to use in gnuplot');
         writeln('defrock data\T4\ 7319 withtimes');
         { writeln('Example 3 : minimal data file with data in volts'); }
         { writeln('defrock data\T4\ 7319 notimes volts'); }
         writeln('');
         writeln('Nothing done this time. Rerun with appropriate parameters.');
     end
     else
     begin
         dataDir := ParamStr(1);
         shotName := ParamStr(2);
         if ParamStr(3) = 'withtimes' then withTimeColumn := true;
         if ParamStr(4) = 'volts' then physicalDataScales := false;

         if dataDir[length(dataDir)] <> '/' then
	 begin
	    dataDir := dataDir + '/';  { Assuming that we are on Linux. }
	 end;
	 writeln('Data directory: '+dataDir);
         writeln('Shot name     : '+shotName);

         startup;

	 destinationdisk[1] := dataDir;
         filename[2] := shotName;
         ChDir(dataDir+'descript');
         readmodel;
         ChDir(homeDir);

         filename[1] := shotName;
         ChDir(dataDir+shotName);
         readdata;
         ChDir(homeDir);

         filename[4] := shotName+'A';
         ChDir(dataDir+shotName);
         writetodisc;
         ChDir(homeDir);

         closedown;
     end;
     writeln('End of defrock run.');
end.

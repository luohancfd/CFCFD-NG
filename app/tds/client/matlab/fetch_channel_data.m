function [t, v] = fetch_channel_data( facilityName, shotName, channelName )
% FETCH_CHANNEL_DATA fetches the time-value pairs.
% Input:
% facilityName : string, e.g. 'T3', 'T4', 'X1', 'X2' or 'X3'
% shotName     : string specifying the base-file-name for the shot.
% channelName  : string specifying the channel number.
% Output:
% Returns the data as two numeric arrays (vectors).

% Peter J. 19-April-03

% Collect the data text as a cell array of strings,
% one per line of the original header.
data_text = fetch_text_from_server( facilityName, shotName, channelName, 'data' );

% Split the strings into time and value components.
count = 0; t = zeros(20000,1); v = zeros(20000,1);
for line_text = data_text
    count = count + 1;
    [time_stamp value] = strread( char(line_text), '%s %s' );
    % At this point, it appears that time_stamp and 
    % value are both cell variables.
    t(count) = sscanf( char(time_stamp), '%f' );
    v(count) = sscanf( char(value), '%f' );
end

% Trim the arrays to the true size of the data.
t = t(1:count);
v = v(1:count);

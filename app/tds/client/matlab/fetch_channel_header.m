function [names, values] = fetch_channel_header( facilityName, shotName, channelName )
% FETCH_CHANNEL_HEADER fetches the channel metadata.
% Input:
% facilityName : string, e.g. 'T3', 'T4', 'X1', 'X2' or 'X3'
% shotName     : string specifying the base-file-name for the shot.
% channelName  : string specifying the channel number.
% Output:
% Returns the name-value pairs as two cell-arrays of strings.

% Peter J. 19-April-03

% Collect the header text as a cell array of strings,
% one per line of the original header.
header_text = fetch_text_from_server( facilityName, shotName, channelName, 'info' );

% Split the strings into name and value components.
count = 0; names = {}; values = {};
for line_text = header_text
    count = count + 1;
    [a, b] = strread( char(line_text), '%s %s');
    names{count} = char(a);
    values{count} = char(b);
end

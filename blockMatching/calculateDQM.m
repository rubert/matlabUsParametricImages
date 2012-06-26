%%%get the RFD file
[fname, path, ~] = uigetfile('*.rfd');

%%%Read all the frames in the file
URI2MATNick(fname, path)

%%%
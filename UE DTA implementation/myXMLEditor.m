function [] = myXMLEditor(xml_file_Java)
%myXMLEditor changes the (.xml)-file due to the additional number of paths.
%[] = myXMLEditor(fileName))
%--------------------------------------------------------------------------
%inputs:
%xml_file_Java: name of the (.xml) network file for Java.
%
% by Manuel Jakob
% 12 April 2012
%==========================================================================

xDoc = xmlread(fullfile(xml_file_Java));
xDocRoot = xDoc.getDocumentElement;

scenario = xDoc.getElementsByTagName('scenario').item(0);
settings = scenario.getElementsByTagName('settings').item(0);
VehicleTypes = settings.getElementsByTagName('VehicleTypes').item(0);

VehicleTypes.appendChild(newLine_element(xDoc));

myXMLwrite(xml_file_Java,xDoc);

end
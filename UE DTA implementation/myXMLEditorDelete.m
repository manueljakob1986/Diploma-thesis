function [] = myXMLEditorDelete(xml_file_Java)
%myXMLEditorDelete deletes all additional number of paths in the
%(.xml)-file.
%[] = myXMLEditorDelete(xml_file_Java)
%--------------------------------------------------------------------------
%inputs:
%xml_file_Java: name of the (.xml) network file for Java.
%
% by Manuel Jakob
% 29 April 2012
%==========================================================================

xDoc = xmlread(fullfile(xml_file_Java));
xDocRoot = xDoc.getDocumentElement;

scenario = xDoc.getElementsByTagName('scenario').item(0);
settings = scenario.getElementsByTagName('settings').item(0);
VehicleTypes = settings.getElementsByTagName('VehicleTypes');

n_types = VehicleTypes.getLength;
for i = 0:n_types-1
    item = VehicleTypes.item(i);
    item.removeChild(item.getLastChild);
end
myXMLwrite(xml_file_Java,xDoc);

end

function out = newLine_element(xDoc)
    nL_element = xDoc.createElement('vehicleType');
    nL_element.setAttribute('name','New Path');
    nL_element.setAttribute('weight','1');
    out = nL_element;
end
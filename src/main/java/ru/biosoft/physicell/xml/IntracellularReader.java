package ru.biosoft.physicell.xml;

import org.w3c.dom.Element;

import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;

public interface IntracellularReader
{
    public void readIntracellular(Element el, Model model, CellDefinition cd)  throws Exception;
}

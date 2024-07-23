package ru.biosoft.physicell.xml;

import java.io.File;

import org.w3c.dom.Element;

import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;

public interface IntracellularReader
{
    public void readIntracellular(File f, Element el, Model model, CellDefinition cd) throws Exception;

    public void readIntracellular(Element el, Model model, CellDefinition cd) throws Exception;
}

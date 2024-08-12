package ru.biosoft.physicell.xml;

import java.io.File;
import java.nio.file.Path;
import java.util.Map;

import org.w3c.dom.Element;

import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;

public interface IntracellularReader
{
    public void setAdditionalFiles(Map<String, File> additional);
    public void readIntracellular(Path path, Element el, Model model, CellDefinition cd) throws Exception;
}

package ru.biosoft.physicell.xml;

import org.w3c.dom.Element;

import ru.biosoft.physicell.core.CellDefinition;

public abstract class FunctionsReader extends ModelReaderSupport
{
    public abstract void readFunctions(Element element, CellDefinition cd);
}

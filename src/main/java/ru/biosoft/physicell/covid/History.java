package ru.biosoft.physicell.covid;

public class History
{
    private final int[] buffer;
    private int head = 0; // points to the oldest element
    private int size = 0;
    private final int capacity;

    public History(int capacity)
    {
        this.capacity = capacity;
        this.buffer = new int[capacity];
    }

    public void addFront(int value)
    {
        head = ( head - 1 + capacity ) % capacity;
        buffer[head] = value;
        if( size < capacity )
        {
            size++;
        }
    }

    public boolean isFull()
    {
        return size == capacity;
    }

    // Get element at back (oldest element)
    public int getBack()
    {
        if( size == 0 )
            throw new IllegalStateException( "Buffer is empty" );
        int backIndex = ( head + size - 1 ) % capacity;
        return buffer[backIndex];
    }
}
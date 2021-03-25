


class initializer
{
    public:
        static initializer& getInstance()
        {
            static initializer    instance; // Guaranteed to be destroyed.
                                  // Instantiated on first use.
            return instance;
        }


    void init()
    {
        if (not initialized)
        {
            amrex::Initialize(MPI_COMM_WORLD);
            initialized=true;
        }
    }

    ~initializer()
    {
        if (initialized)
        {
            amrex::Finalize();
        }
    }

    private:
        initializer() : initialized(false) {}


    public:
        initializer(initializer const&)               = delete;
        void operator=(initializer const&)  = delete;

    bool initialized;

};
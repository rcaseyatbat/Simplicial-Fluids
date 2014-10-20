    std::cout << "Total w : " << totalVorticity << std::endl;
    double fullW = 0;
    for (int i = 0; i < size; i++) {
        fullW += PrimalVertices[i].vorticity;
        if (PrimalVertices[i].boundary == 1) {
            //std::cout << "boundary: " << i << std::endl;
            //std::cout << "v: " << i << " " << PrimalVertices[i].vorticity << std::endl;
        }
    }
    std::cout << "All w:" << fullW << std::endl;

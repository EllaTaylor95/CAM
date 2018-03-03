


n_new = n+1;
        x_new = new double[n_new];

        for (int i = 0; i<k; i++)
            {
            x_new[i] = x[i];
            std::cout<<"x_new["<<i<<"] = "<<x_new[i]<<"\n";
            }

        x_new[k] = 0.5*(x[k-1]+x[k]);
        std::cout<<"x_new["<<k<<"] = "<<x_new[k]<<"\n";

        for (int i = k; i<n; i++)
            {
            x_new[i+1] = x[i];
            std::cout<<"x_new["<<i+1<<"] = "<<x_new[i]<<"\n";
            }
        }
    else
    {
    n_new = n;
    std::cout<<"n_new = "<<n<<"\n";
    x_new = new double[n_new];

    for (int i = 0; i<n_new; i++)
        {
        x_new[i] = x[i];
        std::cout<<"x_new["<<i<<"] = "<<x_new[i]<<"\n";
        }
    }
    n=n_new;
    delete[] x;
    x = new double[n_new];
    for(int i =0; i<n_new; i++)
        {
        x[i] = x_new[i];
        std::cout<<"x["<<i<<"] = "<<x[i]<<"\n";
        }
    delete[] x_new;

    std::cout<< counter <<std::setw(10)<< n-1 <<std::setw(20)<<Gauss_quad1<<std::setw(30)<<Gauss_quad3<<std::setw(40)<<Est<<std::endl;
    }

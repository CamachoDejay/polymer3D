# Fast and Precise 3D Fluorophore Localization based on Gradient Fitting
### 2015 Scientific reports
Ma, H., Xu, J., Jin, J., Gao, Y., Lan, L., & Liu, Y. Scientific Reports, 5(1), 14335. [link](https://doi.org/10.1038/srep14335)




$$
D = \sum_{m,n} \frac{(e(x_c-m) g_y-(y_c-n) g_x )^2  }{(e_0^2 (x_0-m)^2+(y_0-n)^2 )(g_x^2+g_y^2 ) }W
$$

where, $e_0$ and $(x_0,y_0)$ are the initial values of $e$ and $(x_c,y_c)$ estimated by the centroid method (see [QuickPALM](link)).:
$$
W= (g_x^2+g_y^2 )* \sqrt{(x_0-m)^2+(y_0-n)^2}
$$

Now the task is to solve the system of equations:
$$
\frac{\partial D}{\partial x_c}=0 \\
\frac{\partial D}{\partial y_c}=0 \\
\frac{\partial D}{\partial e}=0


$$

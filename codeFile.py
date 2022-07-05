import numpy as np
import cmath as cm
import matplotlib.pyplot as plt

# Writing a function that uses convolution(by mirror padding) to denoise a signal using the kernel h
def denoise(y):
    h = np.full((1,5),(1/5))
    denoised_signal = []
    for i in range(len(y)):
        averaged_sample = 0
        for k in range(-2,3):
            if i < 2:
                if i + k >= 0:
                    averaged_sample += y[i + k] / 5
                else:
                    averaged_sample += y[i - k] / 5
            if 2 <= i <= len(y) - 3:
                averaged_sample += y[i + k] / 5
            if i > len(y) - 3:
                if i + k <= len(y) - 1:
                    averaged_sample += y[i + k] / 5
                else:
                    averaged_sample += y[i - k] / 5
        denoised_signal.append(round(averaged_sample.real,4)+round(averaged_sample.imag,4)*1j)
    return denoised_signal


# A function that calculates the DTFT of a sampled signal and returns 193 samples
def DFT(y):
    d = []
    N = 193
    for k in range(N):
        samples_of_dft = 0
        for n in range(len(y)):
            samples_of_dft+=y[n]*(cm.exp((-1j)*(2*cm.pi/N)*k*n))
        d.append(samples_of_dft)
    return d


# A function that calculates the descrete inverse fourior transform of a signal and returns 193 samples
def INV_DFT(z):
    recovered_signal = []
    number_of_samples_of_inverse_DTFT = 193
    for k in range(number_of_samples_of_inverse_DTFT):
        sample_value = 0
        for n in range(len(z)):
            sample_value+=(z[n]*(cm.exp(2*cm.pi*1j*k*n/number_of_samples_of_inverse_DTFT)))/number_of_samples_of_inverse_DTFT
        recovered_signal.append(round(sample_value.real,4)+round(sample_value.imag,4)*1j)
    return recovered_signal


# h is the blurring kernel
h=[1/16,4/16,6/16,4/16,1/16]


# as we know that h[0]=6/16 so in order to compensate the timeshift we need to use the timeshift property of DTFT
original_DTFT_of_h=[]
for i in range(len(DFT(h))):
    z = cm.exp((4*cm.pi/193)*(1j)*i)*DFT(h)[i]
    original_DTFT_of_h.append(round(z.real,4)+round(z.imag,4)*1j)
# now as we can see the DTFT consists of many terms that tend to 0
# so in order to avoid devision by 0 we can actually clip the lower bound of the original_DTFT_of_h
# by hit and trial we can find that if we clip the lower bound to 0.35 then it is sufficiant to deblur the signal.
# hence the code for the same is as follows
for i in range(len(original_DTFT_of_h)):
    if original_DTFT_of_h[i].real <= 0.35:
        original_DTFT_of_h[i] = 0.35 + 0j



# application time!!!  PHEW!!!
print("INPUT THE SAMPLED SIGNAL Y:")
y = []
for i in range(193):
    y.append(float(input()))

# recovering signal by first denoising and then deblurring
denoised_y = denoise(y)
a = [] # array of samples of DFT(denoised_y)/original_DTFT_of_h
for i in range(193):
    c = DFT(denoised_y)[i] / original_DTFT_of_h[i]
    a.append(round(c.real, 4) + round(c.imag, 4) * 1j)
x_1 = INV_DFT(a) # recovered signal x_1[n]
print("THE SIGNAL x_1[n] is:")
print(np.array(x_1))

#recovering signal by first deblurring and then denoising
b=[]
for i in range(193):
    c = DFT(y)[i] / original_DTFT_of_h[i]
    b.append(round(c.real, 4) + round(c.imag, 4) * 1j)
inv_f = INV_DFT(b)
x_2 = denoise(inv_f)
print("""

THE SIGNAL x_2[n] is:""")
print(np.array(x_2))


print("""
INPUT THE SAMPLED X SIGNAL(ORIGINAL SIGNAL): """)
x = []
for i in range(193):
    x.append(float(input()))


d = [] #array containig absolute error in samples of x_1 and X(original signal)
for i in range(193):
    d.append(abs(x_1[i] - x[i]))
print("the average error in the signal x_1 is:",end=' ')
print(str(sum(np.array(d))/193))


e = [] #array containig absolute error in samples of x_2 and X(original signal)
for i in range(193):
    e.append(abs(x_2[i]-x[i]))
print("the average error in the signal x_2 is:",end=' ')
print(str(sum(np.array(e))/193))

# an array of time
t = np.arange(193)

#plotting x with time
plt.plot(t,x)
plt.xlabel("t")
plt.ylabel("x")
plt.show()

#plotting x_1 with time
plt.plot(t,x_1)
plt.xlabel("t")
plt.ylabel("x_1")
plt.show()

#plotting x_2 with time
plt.plot(t,x_2)
plt.xlabel("t")
plt.ylabel("x_2")
plt.show()
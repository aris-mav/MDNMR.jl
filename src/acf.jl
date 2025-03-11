export ACF
"""
Autocorrelation function
 FFT based recipe:
1.	Pad vector-a by an equal number of zeros. Thus, [1 2 3 0 0 0]
2.	Take the discrete FFT of the new array. Call this F.
3.	Take the conjugate. Call this F*
4.	Compose F times F* (you should do term by term multiplication). Call this Cff.
5.	Take the inverse FFT of Cff, and take only the first 3 terms. Call this ACF.
6.	Now normalize by the vector [3, 2, 1]. That is your answer.
"""
function ACF(v::AbstractVector)

    a = [v ; zeros(length(v))]
    F = fft(a)
    Cff = F .* vec(F')
    ACF = real.(ifft(Cff)[1:length(v)] ./ reverse(collect(1:length(v))))

end


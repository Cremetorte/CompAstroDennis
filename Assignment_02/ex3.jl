using PrettyTables

function midnightformula(a, b, c)
    sol_1 = (-b + sqrt(b^2 - 4*a*c))/2/a  # "normal" approach
    sol_2 = (-b - sqrt(b^2 - 4*a*c))/2/a
    return (sol_1, sol_2)
end


function alternative_midnight(a,b,c)
    sol_1 = 2*c / (-b - sqrt(b^2-4a*c))  # as in the PDF
    sol_2 = (-b - sqrt(b^2 - 4*a*c))/2/a     # as before
    return (sol_1, sol_2)
end


function automatic_midnight(a,b,c)
    if b == sqrt(b^2 - 4*a*c)
        return alternative_midnight(a,b,c)
    else
        return midnightformula(a,b,c)
    end
end




c = [10.0^(-i) for i in 14:20]
a, b = 1.0, 1.0




# Normal mode with catastrophic cancellation
solutions = midnightformula.(a, b, c)
println("Lösungen mit gewöhnlicher Mitternachtsformel")
display(solutions)
#=
 (-9.992007221626409e-15, -0.99999999999999)
 (-9.992007221626409e-16, -0.999999999999999)
 (-1.1102230246251565e-16, -0.9999999999999999) <-- Maschinenpräzision liegt in dieser Größenordnung
 (0.0, -1.0)                                    <-- deshalb wird x_1 hier 0 und x_2 wird 1.0
 (0.0, -1.0)                                    <-- bleibt ab hier
 (0.0, -1.0)
 (0.0, -1.0)
=#


# Adapted mode without catastrophic cancellation
solutions = alternative_midnight.(a, b, c)
println("\nLösungen mit angepasster Mitternachtsformel")
display(solutions)
#=
 (-1.00000000000001e-14, -0.99999999999999)
 (-1.000000000000001e-15, -0.999999999999999)
 (-1.0000000000000001e-16, -0.9999999999999999) <-- Ungenauere Ergebnisse bis zur Maschinenpräzision
 (-1.0e-17, -1.0)                               <-- Hier bleiben jetzt aber die Größenordnungen erhalten
 (-1.0e-18, -1.0)
 (-1.0e-19, -1.0)
 (-1.0e-20, -1.0)
=#


# auto mode: checks if b = √(b^2-4ac)
solutions = automatic_midnight.(a, b, c)
println("\nLösungen mit automatischer Mitternachtsformel")
display(solutions)
#=
 (-9.992007221626409e-15, -0.99999999999999)
 (-9.992007221626409e-16, -0.999999999999999)
 (-1.1102230246251565e-16, -0.9999999999999999)  <-- Bis hier normal
 (-1.0e-17, -1.0)                                <-- ab hier angepasste Methode
 (-1.0e-18, -1.0)
 (-1.0e-19, -1.0)
 (-1.0e-20, -1.0)
=#
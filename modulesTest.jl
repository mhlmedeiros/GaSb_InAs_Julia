module PrintHello
    export foo, spam
    foo() = println("Hello")
    foo(name :: String) = print("Hello ", name)
    spam(n :: Int64) = print("Spam"^n)
end  # module PrintHello

module PrintHi
    export greeting
    greeting() = println("Hi")
end  # module PrintHi

# import .PrintHello
#
# using .PrintHello
#
# foo()

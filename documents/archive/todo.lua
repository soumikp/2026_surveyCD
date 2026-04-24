local todos = {}

-- Pass 1: Collect the notes and format them inline
local collect_todos = {
  Span = function(el)
    if el.classes:includes('todo') then
      -- Collect the text string of the to-do
      table.insert(todos, pandoc.utils.stringify(el))

      -- Apply red coloring for HTML outputs
      if FORMAT == "html" then
        el.attributes['style'] = 'color: red; font-weight: bold;'
        return el
        
      -- Apply red coloring for PDF/LaTeX outputs
      elseif FORMAT:match('latex') or FORMAT:match('pdf') then
        local colored_content = {pandoc.RawInline('latex', '\\textcolor{red}{\\textbf{')}
        for _, v in ipairs(el.content) do 
          table.insert(colored_content, v) 
        end
        table.insert(colored_content, pandoc.RawInline('latex', '}}'))
        return pandoc.Span(colored_content)
      end
    end
  end
}

-- Pass 2: Inject the compiled list at your placeholder
local inject_todos = {
  Div = function(el)
    if el.identifier == 'list-of-todos' then
      -- If there are no todos, return nothing
      if #todos == 0 then 
          return pandoc.Null() 
      end

      -- Build the bulleted list
      local list_items = {}
      for _, t in ipairs(todos) do
        table.insert(list_items, {pandoc.Plain({pandoc.Str(t)})})
      end

      -- Return a header and the list
      return pandoc.Div({
        pandoc.Header(3, {pandoc.Str("To-Do List")}, pandoc.Attr("", {"unnumbered"})),
        pandoc.BulletList(list_items)
      })
    end
  end
}

-- Return the passes in the order they should be executed
return {collect_todos, inject_todos}